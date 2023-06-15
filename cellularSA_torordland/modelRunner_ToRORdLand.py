import numpy as np
from model_ToRORd_Land import model_ToRORd_Land
def modelRunner_ToRORdLand(X0, options, parameters, beats, ignoreFirst):


    cellType = 0
    if (hasattr(parameters, 'cellType')):
        cellType = parameters.cellType


    # if the number of simulated beat is to be printed out.
    verbose = False
    if (hasattr(parameters, 'verbose')):
        verbose = parameters.verbose
    

    nao = 140
    if (hasattr(parameters, 'nao')):
        nao = parameters.nao
    
    cao = 1.8
    if (hasattr(parameters, 'cao')):
        cao = parameters.cao
    
    ko = 5
    if (hasattr(parameters, 'ko')):
        ko = parameters.ko
    
    ICaL_fractionSS = 0.8
    if (hasattr(parameters, 'ICaL_fractionSS')):
        ICaL_fractionSS = parameters.ICaL_fractionSS
    
    INaCa_fractionSS = 0.35
    if (hasattr(parameters, 'INaCa_fractionSS')):
        INaCa_fractionSS = parameters.INaCa_fractionSS
    

    INa_Multiplier = 1
    if (hasattr(parameters, 'INa_Multiplier')):
        INa_Multiplier = parameters.INa_Multiplier
    
    ICaL_Multiplier = 1
    if (hasattr(parameters, 'ICaL_Multiplier')):
        ICaL_Multiplier = parameters.ICaL_Multiplier
    
    Ito_Multiplier = 1
    if (hasattr(parameters, 'Ito_Multiplier')):
        Ito_Multiplier = parameters.Ito_Multiplier
    
    INaL_Multiplier = 1
    if (hasattr(parameters, 'INaL_Multiplier')):
        INaL_Multiplier = parameters.INaL_Multiplier
    
    IKr_Multiplier = 1
    if (hasattr(parameters, 'IKr_Multiplier')):
        IKr_Multiplier = parameters.IKr_Multiplier
    
    IKs_Multiplier = 1
    if (hasattr(parameters, 'IKs_Multiplier')):
        IKs_Multiplier = parameters.IKs_Multiplier
    
    IK1_Multiplier = 1
    if (hasattr(parameters, 'IK1_Multiplier')):
        IK1_Multiplier = parameters.IK1_Multiplier
    
    IKCa_Multiplier = 0
    if (hasattr(parameters, 'IKCa_Multiplier')):
        IKCa_Multiplier = parameters.IKCa_Multiplier
    
    IKb_Multiplier = 1
    if (hasattr(parameters, 'IKb_Multiplier')):
        IKb_Multiplier = parameters.IKb_Multiplier
    
    INaCa_Multiplier = 1
    if (hasattr(parameters, 'INaCa_Multiplier')):
        INaCa_Multiplier = parameters.INaCa_Multiplier
    
    INaK_Multiplier = 1
    if (hasattr(parameters, 'INaK_Multiplier')):
        INaK_Multiplier = parameters.INaK_Multiplier
    
    INab_Multiplier = 1
    if (hasattr(parameters, 'INab_Multiplier')):
        INab_Multiplier = parameters.INab_Multiplier
    
    ICab_Multiplier = 1
    if (hasattr(parameters, 'ICab_Multiplier')):
        ICab_Multiplier = parameters.ICab_Multiplier
    
    IpCa_Multiplier = 1
    if (hasattr(parameters, 'IpCa_Multiplier')):
        IpCa_Multiplier = parameters.IpCa_Multiplier
    
    ICaCl_Multiplier = 1
    if (hasattr(parameters, 'ICaCl_Multiplier')):
        ICaCl_Multiplier = parameters.ICaCl_Multiplier
    
    IClb_Multiplier = 1
    if (hasattr(parameters, 'IClb_Multiplier')):
        IClb_Multiplier = parameters.IClb_Multiplier
    
    Jrel_Multiplier = 1
    if (hasattr(parameters, 'Jrel_Multiplier')):
        Jrel_Multiplier = parameters.Jrel_Multiplier
    
    Jup_Multiplier = 1
    if (hasattr(parameters, 'Jup_Multiplier')):
        Jup_Multiplier = parameters.Jup_Multiplier
    
    aCaMK_Multiplier = 1
    if (hasattr(parameters, 'aCaMK_Multiplier')):
        aCaMK_Multiplier = parameters.aCaMK_Multiplier
    
    taurelp_Multiplier = 1
    if (hasattr(parameters, 'taurelp_Multiplier')):
        taurelp_Multiplier = parameters.taurelp_Multiplier
    
    kws_Multiplier = 1
    if (hasattr(parameters, 'kws_Multiplier')):
        kws_Multiplier = parameters.kws_Multiplier
    
    kuw_Multiplier = 1
    if (hasattr(parameters, 'kuw_Multiplier')):
        kuw_Multiplier = parameters.kuw_Multiplier
    
    ksu_Multiplier = 1
    if (hasattr(parameters, 'ksu_Multiplier')):
        ksu_Multiplier = parameters.ksu_Multiplier
    
    ca50_Multiplier = 1
    if (hasattr(parameters, 'ca50_Multiplier')):
        ca50_Multiplier = parameters.ca50_Multiplier
    

    extraParams = []
    if (hasattr(parameters, 'extraParams')):
        extraParams = parameters.extraParams
    

    vcParameters = []
    if (hasattr(parameters, 'vcParameters')):
        vcParameters = parameters.vcParameters
    
    apClamp = []
    if (hasattr(parameters, 'apClamp')):
        apClamp = parameters.apClamp
    

    stimAmp = -53
    if (hasattr(parameters, 'stimAmp')):
        stimAmp = parameters.stimAmp
    
    stimDur = 1
    if (hasattr(parameters, 'stimDur')):
        stimDur = parameters.stimDur
    # Ca50 = 0.805 if (hasattr(parameters, 'Ca50')):
    Ca50 = parameters.Ca50
    
    trpnmax = 0.07
    if (hasattr(parameters, 'trpnmax')):
        trpnmax = parameters.trpnmax

    CL = parameters.bcl
    time = cell(beats, 1)
    X = cell(beats, 1)
    for n in range(beats):
        if (verbose):
            print('Beat = ', str(n))
    
        if ((exist('ode15sTimed')) == 2): # if timed version provided, it is preferred
        [time{n}, X{n}] = ode15sTimed(parameters.model, [0 CL], X0, options, 1, cellType, ICaL_Multiplier, ...
        INa_Multiplier, Ito_Multiplier, INaL_Multiplier, IKr_Multiplier, IKs_Multiplier, IK1_Multiplier, IKCa_Multiplier, IKb_Multiplier, INaCa_Multiplier, ...
        INaK_Multiplier, INab_Multiplier, ICab_Multiplier, IpCa_Multiplier, ICaCl_Multiplier, IClb_Multiplier, Jrel_Multiplier, Jup_Multiplier, aCaMK_Multiplier, ...
        taurelp_Multiplier, kws_Multiplier, kuw_Multiplier, ksu_Multiplier, ca50_Multiplier, nao, cao, ko, ICaL_fractionSS, INaCa_fractionSS, stimAmp, stimDur, trpnmax, vcParameters, apClamp, extraParams)
        else
        [time{n}, X{n}] = ode15s(parameters.model, [0 CL], X0, options, 1, cellType, ICaL_Multiplier, ...
        INa_Multiplier, Ito_Multiplier, INaL_Multiplier, IKr_Multiplier, IKs_Multiplier, IK1_Multiplier, IKCa_Multiplier, IKb_Multiplier, INaCa_Multiplier, ...
        INaK_Multiplier, INab_Multiplier, ICab_Multiplier, IpCa_Multiplier, ICaCl_Multiplier, IClb_Multiplier, Jrel_Multiplier, Jup_Multiplier, aCaMK_Multiplier, ...
        taurelp_Multiplier, kws_Multiplier, kuw_Multiplier, ksu_Multiplier, ca50_Multiplier, nao, cao, ko, ICaL_fractionSS, INaCa_fractionSS, stimAmp, stimDur, trpnmax, vcParameters, apClamp, extraParams)


    if isequal(time{n}, -1) # unfinished (unstable) computation - we  here.
    try
        time(1: ignoreFirst) = []
        X(1: ignoreFirst) = []
    catch
    time = []
    X = []
    
    parameters.isFailed = 1
    return
    

    X0 = X
    {n}(size(X
    {n}, 1),:)

    time(1: ignoreFirst) = []
    X(1: ignoreFirst) = []
    return time, X, parameters

