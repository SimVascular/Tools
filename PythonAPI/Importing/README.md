 # Tools for Installing, Uninstalling, and Checking for SV Python Packages
Currently this code serves as a template to be copied into scripts that are
loaded into the SV Python API Interface when additional Python modules are
needed but not included within the base SV Python Interpreter. For reference,
the python modules included with the base SV Python Interpreter can be listed
from the following command.

```python
help("modules")
```

Albeit after much complaining, the SV Python Interpreter will return a
comprehensive list of all the included modules which are available by
default. If your script requires a new or different version of one of
these modules, *be careful*. Changing one of these modules may affect
performance of internal routines and cause unpredicted behaviors to
occur. This is especially true for the `vtk` module which is called
by `sv` specific modules to return PolyData objects. This list of included
modules is also provided here for reference.

|      |      |     |     |
|:----:|:----:|:---:|:---:|
|Cython             | codeop             | netrc              | sys|
|IPython            |collections         | nntplib            | sysconfig|
|__future__         | colorama           | notebook           | tabnanny|
|_ast               | colorsys           | nt                 | tarfile|
|_bisect            | commctrl           | ntpath             | telnetlib|
|_bootlocale        | compileall         | ntsecuritycon      | tempfile|
|_cffi_backend      | concurrent         | nturl2path         | tensorboard|
|_codecs            | configparser       | numbers            | tensorflow|
|_codecs_cn         | contextlib         | numpy              | tensorflow_estimator|
|_codecs_hk         | copy               | odbc               | termcolor|
|_codecs_iso2022    | copyreg            | opcode             | terminado|
|_codecs_jp         | crypt              | operator           | testpath|
|_codecs_kr         | csv                | optparse           | tests|
|_codecs_tw         | ctypes             | os                 | textwrap|
|_collections       | curses             | packaging          | this|
|_collections_abc   | cycler             | pandocfilters      | threading|
|_compat_pickle     | cygwin_bash_kernel | parser             | time|
|_compression       | cython             | parso              | timeit|
|_csv               | cythonmagic        | pasta              | timer|
|_ctypes            | datetime           | pathlib            | tkinter|
|_ctypes_test       | dateutil           | pdb                | token|
|_datetime          | dbi                | perfmon            | tokenize|
|_decimal           | dbm                | pickle             | tornado|
|_distutils_hack    | dde                | pickleshare        | trace|
|_dummy_thread      | decimal            | pickletools        | traceback|
|_elementtree       | decorator          | pip                | tracemalloc|
|_functools         | defusedxml         | pipes              | traitlets|
|_hashlib           | difflib            | pkg_resources      | tty|
|_heapq             | dis                | pkgutil            | turtle|
|_imp               | distutils          | platform           | turtledemo|
|_io                | doctest            | plistlib           | types|
|_json              | dummy_threading    | poplib             | typing|
|_locale            | easy_install       | posixpath          | unicodedata|
|_lsprof            | echo_kernel        | powershell_kernel  | unittest|
|_markupbase        | email              | pprint             | urllib|
|_md5               | encodings          | profile            | urllib3|
|_msi               | ensurepip          | prometheus_client  | uu|
|_multibytecodec    | entrypoints        | prompt_toolkit     | uuid|
|_multiprocessing   | enum               | pstats             | venv|
|_opcode            | errno              | pty                | vtk|
|_operator          | faulthandler       | pvectorc           | warnings|
|_osx_support       | filecmp            | pwlf               | wave|
|_overlapped        | fileinput          | pyDOE              | wcwidth|
|_pickle            | fnmatch            | py_compile         | weakref|
|_pydecimal         | formatter          | pyclbr             | webbrowser|
|_pyio              | fractions          | pycparser          | webencodings|
|_pyrsistent_version| ftplib             | pydoc              | werkzeug|
|_random            | functools          | pydoc_data         | wheel|
|_sha1              | gast               | pyexpat            | widgetsnbextension|
|_sha256            | gc                 | pygments           | win2kras|
|_sha512            | genericpath        | pyparsing          | win32api|
|_signal            | getopt             | pyrsistent         | win32clipboard|
|_sitebuiltins      | getpass            | pythoncom          | win32com|
|_socket            | gettext            | pywin              | win32con|
|_sqlite3           | glob               | pywin32_bootstrap  | win32console|
|_sre               | grpc               | pywin32_testutil   | win32cred|
|_ssl               | gzip               | pywintypes         | win32crypt|
|_stat              | h5py               | pyximport          | win32cryptcon|
|_string            | hashlib            | qtconsole          | win32event|
|_strptime          | heapq              | qtpy               | win32evtlog|
|_struct            | hmac               | queue              | win32evtlogutil|
|_symtable          | html               | quopri             | win32file|
|_testbuffer        | http               | random             | win32gui|
|_testcapi          | idlelib            | rasutil            | win32gui_struct|
|_testimportmultiple| idna               | re                 | win32help|
|_testmultiphase    | imaplib            | regcheck           | win32inet|
|_thread            | imghdr             | regutil            | win32inetcon|
|_threading_local   | imp                | reprlib            | win32job|
|_tkinter           | importlib          | requests           | win32lz|
|_tracemalloc       | importlib_metadata | rlcompleter        | win32net|
|_warnings          | inspect            | rmagic             | win32netcon|
|_weakref           | io                 | run                | win32pdh|
|_weakrefset        | ipaddress          | runpy              | win32pdhquery|
|_win32sysloader    | ipykernel          | sched              | win32pdhutil|
|_winapi            | ipykernel_launcher | scipy              | win32pipe|
|_winxptheme        | ipython_genutils   | select             | win32print|
|_yaml              | ipywidgets         | selectors          | win32process|
|abc                | isapi              | send2trash         | win32profile|
|absl               | itertools          | servicemanager     | win32ras|
|adodbapi           | jedi               | setuptools         | win32rcparser|
|afxres             | jinja2             | shelve             | win32security|
|aifc               | json               | shlex              | win32service|
|antigravity        | jsonschema         | shutil             | win32serviceutil|
|argon2             | jupyter            | signal             | win32timezone|
|argparse           | jupyter_client     | simvascular_python_kernel | win32trace|
|array              | jupyter_console    | simvascular_tcl_kernel | win32traceutil|
|ast                | jupyter_core       | site               | win32transaction|
|astor              | keras_applications | six                | win32ts|
|asynchat           | keras_preprocessing| smtpd              | win32ui|
|asyncio            | keyword            | smtplib            | win32uiole|
|asyncore           | kiwisolver         | sndhdr             | win32verstamp|
|atexit             | lib2to3            | socket             | win32wnet|
|attr               | linecache          | socketserver       | win_unicode_console|
|audioop            | locale             | sqlite3            | winerror|
|autoreload         | logging            | sre_compile        | winioctlcon|
|backcall           | lzma               | sre_constants      | winnt|
|base64             | macpath            | sre_parse          | winperf|
|bdb                | macurl2path        | ssl                | winpty|
|binascii           | mailbox            | sspi               | winreg|
|binhex             | mailcap            | sspicon            | winxpgui|
|bisect             | markdown           | stat               | winxptheme|
|bleach             | markupsafe         | statistics         | wrapt|
|builtins           | marshal            | storemagic         | wsgiref|
|bz2                | math               | string             | xdrlib|
|cProfile           | mimetypes          | stringprep         | xml|
|calendar           | mistune            | struct             | xmlrpc|
|certifi            | mmap               | subprocess         | xxsubtype|
|cffi               | mmapfile           | sunau              | yaml|
|cgi                | mmsystem           | sv                 | zipapp|
|cgitb              | modulefinder       | sv_ml              | zipfile|
|chardet            | msilib             | sv_rom_extract_results | zipimport|
|chunk              | msvcrt             | sv_rom_simulation  | zipp|
|cmath              | multiprocessing    | sv_vis             | zlib|
|cmd                | nbconvert          | symbol             | zmq|
|code               | nbformat           | sympyprinting      |    |
|codecs             | netbios            | symtable           |    |

Again, this is meant to be a convenient list of modules that are most likely
available with the base SV Python Interpreter; however, this should be confirmed
on your own machine before assuming that is or is not available.

To check the version of a given module the following code can be run. This may
be important when checking for potential conflict between SV and external
package dependencies. If there is a conflict between versions, always prioritize
the SV package version dependency.

```python
import pkg_resources
pkg_resources.get_distribution('<package_name>').version
```

Example:
```python
>>> import pkg_resources
>>> pkg_resources.get_distribution('six').version
'1.15.0'
```

This tool folder contains the `install`, `uninstall`, and `search` functions.
These functions attempt to find packages and add or remove them to the
SV Python 3.5.5 Interpreter.

## Help

For more information on function usage, inputs, and return values, run the
`help(<function_name>)` command on a given function.

```python
help(install)
help(uninstall)
help(search)
```
