from pip._internal import main as pipmain

def install(package_name,version=None,inequality=None):
    """
    Install a python package from within the python SV Interpreter.

    Parameters
    ----------

    package_name : str
                   python package name for lookup on https://pypi.org/
    version : str, optional
              package version to install. Default is None meaning that
              no constraints on package version will be supplied
    inequality : str, optional
              supplied as an additional argument if multiple packages
              satisfy the desired requirements above or below the given
              package version specified. By default inequality is None,
              meaning that the inequality argument is ignored; however
              if a version is supplied without an accompanying inequality
              argument then inequality defaults to strict equality '=='

              Accepted string arguments are the following:
              '==', '<', '>', '<=', '>='

    Returns
    -------

    success : bool
              if a package was successfully installed success will be True,
              otherwise success will return False
    """
    # Sanitizing Inputs for package requirements
    if not isinstance(package_name,str):
        print('package_name argument must be of type -> str')
        return False
    if not isinstance(version,str) and not version is None:
        print('version argument must be of type -> str')
        return False
    if not isinstance(inequality,str) and inequality is not None:
        print('inequality argument must be of type -> str')
        return False
    # Checking logical package constraints (if applicable)
    if inequality is not None:
        if inequality is not in ['==','<','>','<=','>=']:
            print('inequality must be one of the following: == < > <= >=')
            return False
        if version is None:
            print('version must be given if an inequality argument is given')
            return False
    # Attempting package fetching and installation
    try:
        if version is None:
            pipmain(['install',package_name])
        elif version is not None and inequality is None:
            pipmain(['install',package_name+'=='+version])
        elif version is not None and inequality is not None:
            pipmain(['install',package_name+inequality+version])
        success = True
    except:
        print('package installation failed')
        success = False
    return success
