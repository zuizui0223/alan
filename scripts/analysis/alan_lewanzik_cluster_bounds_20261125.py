"""Published aggregate sensitivity, Lewanzik & Voigt 2014.
14 plants, 63 marked fruits per arm, dark 63 removed, lit 49 removed;
lit plants have 1--10 marked fruits, reported min and max attained.
Enumeration does not reconstruct actual plant-level observations.
"""
from fractions import Fraction
from functools import lru_cache
@lru_cache(None)
def solve(k,n,f,mn,mx):
    if k==0:return (Fraction(0),Fraction(0)) if n==f==0 and mn and mx else None
    if not k<=n<=10*k or not 0<=f<=n:return None
    lo=hi=None
    for a in range(1,11):
        for b in range(min(a,f)+1):
            z=solve(k-1,n-a,f-b,mn or a==1,mx or a==10)
            if z is None:continue
            l,h=z;v=Fraction(b,a);l+=v;h+=v
            lo=l if lo is None else min(lo,l)
            hi=h if hi is None else max(hi,h)
    return (lo,hi) if lo is not None else None
if __name__=="__main__":
    z=solve(14,63,14,False,False)
    assert z==(Fraction(7,5),Fraction(91,10))
    print("Plant-equal removal proportion identification bounds:",float(1-z[1]/14),float(1-z[0]/14))
    # Reported paired t, mean difference and n=12; rounded published statistics.
    se=112/4.1;t=2.200985160082949
    print("Conditional delay difference, SE and approximate t11 interval:",112,se,[112-t*se,112+t*se])
