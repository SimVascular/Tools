#!/usr/bin/env python

# author: Zachary Sexton
#   date: 8-31-21

from pip._internal import main as pipmain

def uninstall(package_name):
    """
    Install a python package from within the python SV Interpreter.

    Parameters
    ----------

    package_name : str
                   python package name for lookup within the site-packages
                   directory for the SV Interpreter

    Returns
    -------

    success : bool
              if a package was successfully uninstalled success will be True,
              otherwise success will return False
    """
    # Sanitizing Inputs for package requirements
    if not isinstance(package_name,str):
        print('package_name argument must be of type -> str')
        return False
    try:
        pipmain(['uninstall',package_name])
        success = True
    except:
        print('package uninstallation failed')
        success = False
    return success
