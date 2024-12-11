class HealthyBiomarkerRanges:
    def __init__(self):
        self.healthy_ranges = {'EDVL': [88, 161],
                                # [mL] female UKB https://jcmr-online.biomedcentral.com/articles/10.1186/s12968-017-0327-9/tables/3 95% confidence interval
                                'ESVL': [31, 68],
                                # [mL] female UKB https://jcmr-online.biomedcentral.com/articles/10.1186/s12968-017-0327-9/tables/3
                                'SVL': [49, 100],
                                # [mL] female UKB https://jcmr-online.biomedcentral.com/articles/10.1186/s12968-017-0327-9/tables/3
                                'LVEF': [51, 70],
                                # % female UKB https://jcmr-online.biomedcentral.com/articles/10.1186/s12968-017-0327-9/tables/3
                                'EDPL': [5000, 15000],
                                # # [Barye] +- standard error of the mean https://link.springer.com/article/10.1007/s12265-018-9816-y/tables/1
                                'PmaxL': [143000, 153000],
                                # [Barye] +- standard error of the mean https://link.springer.com/article/10.1007/s12265-018-9816-y/tables/1
                                'EDVR': [85, 166],
                                # [mL] female UKB https://jcmr-online.biomedcentral.com/articles/10.1186/s12968-017-0327-9/tables/3
                                'ESVR': [27, 77],
                                # [mL] female UKB https://jcmr-online.biomedcentral.com/articles/10.1186/s12968-017-0327-9/tables/3
                                'SVR': [48, 99],
                                # [mL] female UKB https://jcmr-online.biomedcentral.com/articles/10.1186/s12968-017-0327-9/tables/3
                                'RVEF': [47, 68],
                                # % female UKB https://jcmr-online.biomedcentral.com/articles/10.1186/s12968-017-0327-9/tables/3
                                'PmaxR': [28000, 36000],
                                # [Barye] +- standard deviation https://www.sciencedirect.com/science/article/pii/016752739190221A
                                'es_ed_avpd': [1.4, 1.9],
                                # [cm] ES - ED for controls https://journals.physiology.org/doi/full/10.1152/ajpheart.01148.2006 range, lower range comes from
                               # https://heart.bmj.com/content/heartjnl/78/3/230.full.pdf for zero mortality cut off
                                'es_ed_apical_displacement': [-0.001, 0.51],
                                # [cm] ES - ED for controls https://journals.physiology.org/doi/full/10.1152/ajpheart.01148.2006 range
                                'ED_wall_thickness': [0.571, 1.381],
                                # [cm] min and max mean value from 17 AHA regions https://www.sciencedirect.com/science/article/pii/S1361841519300489
                                'ES_wall_thickness': [1.07, 1.749],
                                # [cm] min and max mean value from 17 AHA regions https://www.sciencedirect.com/science/article/pii/S1361841519300489
                                'diff_lv_wall_thickness': [0.084, 0.429],
                                # [cm] min and max mean values from 17 AHA regions https://www.sciencedirect.com/science/article/pii/S1361841519300489
                                'peak_lambda': [1.1, 1.3], # REF NEEDED
                                'min_lambda': [0.85, 0.88], # same as lambda_mid
                                'lambda_endo': [0.85, 0.88],
                                # range https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.28724
                                'lambda_mid': [0.85, 0.88],
                                # range https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.28724
                                'lambda_epi': [0.85, 0.86],
                                # range https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.28724
                                'cross_lambda_endo': [0.79, 0.81],
                                # range https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.28724
                                'cross_lambda_mid': [0.83, 0.86],
                                # range https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.28724
                                'cross_lambda_epi': [0.88, 0.9],
                                # range https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.28724
                                'max_four_chamber_Ell': [-0.16, -0.16],
                                'max_mid_Err': [0.22, 0.29],
                                'max_mid_Ecc': [-0.14, -0.08], # REF NEEDED
                                'qt_dur_mean': [378, 419],
                                # [ms] UKB Imaging + EST combined interquartile range https://www.ahajournals.org/doi/full/10.1161/CIRCGEN.120.003231 (Table 1)
                                'jt_dur_mean' : [294, 335], # [ms] UKB Imaging + EST combined interquartile range https://www.ahajournals.org/doi/full/10.1161/CIRCGEN.120.003231 (Table 1)
                                't_pe_mean': [50, 72],
                                # [ms] LTVA interquartile range UKB https://www.ahajournals.org/action/downloadSupplement?doi=10.1161%2FJAHA.121.025897&file=jah37786-sup-0001-data-tabs-figs.pdf
                                't_dur_mean': [100, 114],
                                # [ms] UKB normal LV mass interquartile range https://academic.oup.com/ehjdh/article/4/4/316/7188159?login=true
                                'qrs_dur_mean': [77, 91],
                                # [ms] UKB Imaging + EST combined interquartile range https://academic.oup.com/ehjdh/article/4/4/316/7188159?login=true
                                'dvdt_ejection': [308, 512], # ml/s, +- std left ventricular dV/dt ejection rate from MRI https://www.sciencedirect.com/science/article/pii/S0894731702744702#fig2
                                'dvdt_filling': [251, 471], # mL/s, +- std left ventricular filling rate from MRI https://www.sciencedirect.com/science/article/pii/S0894731702744702#fig2
                                'hypertension_dvdt_ejection': [222, 430], # mL/s, +- std, ejection rate from echo for HTN patients https://www.sciencedirect.com/science/article/pii/S0894731702744702#fig2
                                'hypertension_dvdt_filling': [109, 327], # mL/s, +- std, filling rate from echo for HTN patients https://www.sciencedirect.com/science/article/pii/S0894731702744702#fig2
                                'dpdt_max_mmHg': [895, 1347], # mmHg/s, +- std, invasively measured in patients after CRT application. https://academic.oup.com/europace/article/13/7/984/446139
                                'dpdt_max': [900000, 2100000], # Barye/s, converted from measurement https://academic.oup.com/europace/article/13/7/984/446139
                                'edpvr_a_klotz': [28.2], # mmHg/mL, https://journals.physiology.org/doi/full/10.1152/ajpheart.01240.2005
                                'edpvr_b_klotz' : [2.79], # dimensionless, https://journals.physiology.org/doi/full/10.1152/ajpheart.01240.2005
                                'edpvr_v_intercept': [0, 0],
                                'espvr': [2.6,3.6] # mmHg/mL, https://www.ahajournals.org/doi/epdf/10.1161/01.CIR.65.5.988
                                }

