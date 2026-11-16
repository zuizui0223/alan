"""Site-level test of modification by baseline open-water cover, with spatial block effects."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.stats import t


def fit(g):
    z=pd.DataFrame({"light":g.light/10,"water":g.area0/10})
    z["interaction"]=z.light*z.water
    controls=["light0","area_trend","freq_trend","temp0","rain0","temp_trend","rain_trend",
              "built0","built","effort","days","spring","summer","autumn","start","span","years","island","nc","nr"]
    outcomes=["candidate_rate","reference_rate","rate_difference"]
    data=pd.concat([z,g[controls+outcomes]],axis=1)
    a=(data-data.groupby(g.block).transform("mean")).to_numpy()
    x,y=a[:,:-3],a[:,-3:]
    b=np.linalg.lstsq(x,y,rcond=1e-10)[0]
    other=np.delete(x,2,axis=1)
    v=x[:,2]-other@np.linalg.lstsq(other,x[:,2],rcond=1e-10)[0]
    ss=v@v
    np.testing.assert_allclose(b[2],v@y/ss)
    score=pd.DataFrame(v[:,None]*(y-x@b)/ss).groupby(g.block.to_numpy()).sum()
    info=pd.Series(v*v/ss).groupby(g.block.to_numpy()).sum()
    k=sum(info>1e-9)
    residual_df=len(g)-len(info)-np.linalg.matrix_rank(x)
    se=np.sqrt((score**2).sum().to_numpy()*k/(k-1)*(len(g)-1)/residual_df)
    half=t.ppf(.975,k-1)*se
    return [dict(outcome=col,sites=len(g),blocks=len(info),effective_blocks=1/(info@info),
                 largest_block_share=info.max(),slope=b[2,j],low=b[2,j]-half[j],high=b[2,j]+half[j])
            for j,col in enumerate(outcomes)]


def main(folder):
    g=pd.read_csv(folder/"coordinate_trends.csv.gz")
    rows=[]
    for label,data in [("all",g),("maintained",g[g.state.eq("maintained")])]:
        rows.extend(dict(scope=label,**row) for row in fit(data))
    out=pd.DataFrame(rows)
    out.to_csv(folder/"open_water_light_interaction.csv",index=False)
    print(out.round(2).to_string(index=False))


if __name__=="__main__":
    main(Path(sys.argv[1]))
