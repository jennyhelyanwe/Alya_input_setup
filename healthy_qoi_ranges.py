class HealthyBiomarkerRanges:
    def __init__(self):
        self.healthy_ranges = {'LVEDV': [88, 161],
                                # [mL] female UKB https://jcmr-online.biomedcentral.com/articles/10.1186/s12968-017-0327-9/tables/3 95% confidence interval
                                'LVESV': [31, 68],
                                # [mL] female UKB https://jcmr-online.biomedcentral.com/articles/10.1186/s12968-017-0327-9/tables/3
                                'LVSV': [49, 100],
                                # [mL] female UKB https://jcmr-online.biomedcentral.com/articles/10.1186/s12968-017-0327-9/tables/3
                                'LVEF': [51, 70],
                                # % female UKB https://jcmr-online.biomedcentral.com/articles/10.1186/s12968-017-0327-9/tables/3
                                'LVEDP': [0.5, 1.5],
                                # # [kPa] +- standard error of the mean https://link.springer.com/article/10.1007/s12265-018-9816-y/tables/1
                                'LVESP': [14.3, 15.3],
                                # [kPa] +- standard error of the mean https://link.springer.com/article/10.1007/s12265-018-9816-y/tables/1
                                'RVEDV': [85, 166],
                                # [mL] female UKB https://jcmr-online.biomedcentral.com/articles/10.1186/s12968-017-0327-9/tables/3
                                'RVESV': [27, 77],
                                # [mL] female UKB https://jcmr-online.biomedcentral.com/articles/10.1186/s12968-017-0327-9/tables/3
                                'RVSV': [48, 99],
                                # [mL] female UKB https://jcmr-online.biomedcentral.com/articles/10.1186/s12968-017-0327-9/tables/3
                                'RVEF': [47, 68],
                                # % female UKB https://jcmr-online.biomedcentral.com/articles/10.1186/s12968-017-0327-9/tables/3
                                'RVESP': [2.8, 3.6],
                                # [kPa] +- standard deviation https://www.sciencedirect.com/science/article/pii/016752739190221A
                                'AVPD': [1.4, 1.9],
                                # [cm] ES - ED for controls https://journals.physiology.org/doi/full/10.1152/ajpheart.01148.2006 range
                                'apical_displacement': [-0.001, 0.51],
                                # [cm] ES - ED for controls ttps://journals.physiology.org/doi/full/10.1152/ajpheart.01148.2006 range
                                'ED_wall_thickness': [0.571, 1.381],
                                # [cm] min and max mean value from 17 AHA regions https://www.sciencedirect.com/science/article/pii/S1361841519300489
                                'ES_wall_thickness': [1.07, 1.749],
                                # [cm] min and max mean value from 17 AHA regions https://www.sciencedirect.com/science/article/pii/S1361841519300489
                                'wall_thickening': [0.084, 0.429],
                                # [cm] min and max mean values from 17 AHA regions https://www.sciencedirect.com/science/article/pii/S1361841519300489
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
                                'QTc': [369.2, 399.2],
                                # [ms] UKB normal LV mass interquartile range https://academic.oup.com/ehjdh/article/4/4/316/7188159?login=true
                                'Tpe': [50, 72],
                                # [ms] LTVA interquartile range UKB https://www.ahajournals.org/action/downloadSupplement?doi=10.1161%2FJAHA.121.025897&file=jah37786-sup-0001-data-tabs-figs.pdf
                                'T_duration': [100, 114],
                                # [ms] UKB normal LV mass interquartile range https://academic.oup.com/ehjdh/article/4/4/316/7188159?login=true
                                'QRS_duration': [81, 97],
                                # [ms] UKB normal LV mass interquartile range https://academic.oup.com/ehjdh/article/4/4/316/7188159?login=true
                                'diastoic_duration': [0.45],
                                'dvdt_ejection': [308, 512], # ml/s, +- std left ventricular dV/dt ejection rate from MRI https://www.sciencedirect.com/science/article/pii/S0894731702744702#fig2
                                'dvdt_filling': [251, 471], # mL/s, +- std left ventricular filling rate from MRI https://www.sciencedirect.com/science/article/pii/S0894731702744702#fig2
                                'hypertension_dvdt_ejection': [222, 430], # mL/s, +- std, ejection rate from echo for HTN patients https://www.sciencedirect.com/science/article/pii/S0894731702744702#fig2
                                'hypertension_dvdt_filling': [109, 327], # mL/s, +- std, filling rate from echo for HTN patients https://www.sciencedirect.com/science/article/pii/S0894731702744702#fig2
                                'dpdt_max_mmHg': [895, 1347], # mmHg/s, +- std, invasively measured in patients after CRT application. https://academic.oup.com/europace/article/13/7/984/446139
                                'dpdt_max': [119, 180] # kPa/s, converted from measurement https://academic.oup.com/europace/article/13/7/984/446139
                               }

