from thermo.interaction_parameters import IPDB
import os

def get_interaction_parameters(CASs):
    r'''
    Reads a modified version of the Dechame/ChemSep database 
    of interaction parameters for the Peng-Robinson equation of state    
    
    Parameters
    ----------
    CASs: list of str
        List of the CAS names for which kij's shall be retrieved

    
    Returns
    -------
    kijs : list of list of float
        Matrix with retrieved BIPs. Default value is 0 if not available. 
    
    Examples
    --------
    >>> CASs = ['7732-18-5', '74-82-8', '124-18-5']
    >>> get_interaction_parameters(CASs)
    [[0.0, 0.5, 0.5], [0.5, 0.0, 0.0411], [0.5, 0.0411, 0.0]]

    References
    ----------
    ...  [1] Amy L. Dill, Harry Kooijman,  DECHEMA Peng-Robinson Parameters, 1996/2009
    ...  [2] Harry Kooijman and Ross Taylor, Interaction Parameters Data (IPD)
    '''
    directory = os.path.dirname(os.path.abspath(__file__))
    filename = 'pr_ram.json'
    path = os.path.join(directory,filename)
    IPDB.load_json(path,'RAM_PR')
    kijs =  IPDB.get_ip_symmetric_matrix(name='RAM_PR', CASs=CASs, ip = 'kij' )
    return kijs

def pseudo_expand_kijs():
    r'''
    Method to expand the kij matrix to include BIPs for pseudo components.
    The following parameters are set as default 

    '''
    pass

if __name__ == '__main__':
    from chemicals import *
    HC_CASs = [ search_chemical('water').CASs,
            search_chemical('methane').CASs,
            search_chemical('ethane').CASs,
            search_chemical('propane').CASs, 
            search_chemical('i-butane').CASs,
            search_chemical('n-butane').CASs,
            search_chemical('i-pentane').CASs,
            search_chemical('n-hexane').CASs,
            search_chemical('heptane').CASs,
            search_chemical('n-octane').CASs,
            search_chemical('n-nonane').CASs,
            search_chemical('n-decane').CASs]

    kij = get_interaction_parameters(HC_CASs)