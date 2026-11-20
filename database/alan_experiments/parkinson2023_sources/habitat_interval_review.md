# Independent derived-interval review

Reviewed report_habitat_intervals.py and independently recalculated from existing trap_units.csv using scalar arithmetic/sample variances (not a fresh original-workbook extraction). All 12 finite intervals and eight zero-mean gates agree. All 144 existing moderation covariance entries agree with independently centred within-cell cross-products (maximum difference 2.22e-16).

For four independent station groups, b=log(L_litt/D_litt)-log(L_rip/D_rip); V=sum(s_i^2/(n_i*m_i^2)); df=V^2/sum(v_i^2/(n_i-1)). The implemented exp(b +/- t_.975,df*sqrt(V)) is the proposed approximate delta/Welch interval. With n=3, this is an exploratory approximation, not an exact small-sample confidence guarantee. No pseudocount or zero-group removal occurs. Shared dark controls contribute covariance across colours; overlapping endpoints also covary. Individual marginal intervals do not require these cross-contrast covariances, but any joint synthesis would.

All-capture ratio-of-ratios (95% approximate interval; df):

| Source colour | Ratio | Lower | Upper | df |
|---|---:|---:|---:|---:|
| green | 2.070097 | 0.451044 | 9.500855 | 4.195354 |
| purple (violet in paper) | 10.748906 | 2.342851 | 49.315550 | 4.872512 |
| red | 0.791367 | 0.168022 | 3.727270 | 5.478811 |
| white | 2.557664 | 0.509192 | 12.847100 | 5.267183 |

## Source/design scope

[Primary Methods 2a–2b](https://pmc.ncbi.nlm.nih.gov/articles/PMC10613537/) were recovered through indexed primary text: one pond, one night (3 June 2022), randomized light assignment, three stations per colour per habitat, approximately 5 m station spacing and 15 m habitat separation. Authors checked illumination between stations below 1 lux. This supports a conditional independent-station calculation, but habitat itself is not independently replicated. Insect movements could still link nearby stations; separation of light exposure does not prove biological independence. Violet was not illuminance-standardized with other colours. Contrasts describe trap capture at littoral versus adjacent riparian placements, not wet versus dry biomes, population decline, direct predation, or separate study votes. Unadjusted multiple intervals should remain exploratory. The cached article.xml is invalid and was not used as full-text provenance.
