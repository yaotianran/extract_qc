import shutil
import hashlib
import os
import io
import os.path as path
from typing import Union
import requests



def download(url, folder = '.', prefix = ''):
    local_filename = path.join(folder, prefix + url.split('/')[-1])
    with requests.get(url, stream=True) as r:
        with open(local_filename, 'wb') as f:
            shutil.copyfileobj(r.raw, f)

    return local_filename

def md5(file: Union[str, io.BufferedReader]) -> str:
    '''
    A simple hashlib.md5 wrapper. Calculate md5 hash of a file.

    Do remember that the security of the MD5 has been severely compromised since 2012, so it should NOT be used for security purposes.

    Parameters:
        **file**: String.
            The file you want to do the hash

    Return:
        **string**:a 32 bytes MD5 hash string.
    '''


    if not isinstance(file, str) and not isinstance(file, io.BufferedReader):
        message = 'Parameter "file" needs to be a string or a file object, instead a {} is given.'.format(type(file))
        raise TypeError(message)

    if isinstance(file, str):
        file_str = path.realpath(path.expanduser(file))
        if path.isdir(file_str):
            message = 'Parameter "file" needs to be a file, but {} is folder.'.format(file_str)
            raise IsADirectoryError(message)

        if not os.access(file_str, os.R_OK):
            message = "File {} doesn't exists or is not readable.".format(file_str)
            raise PermissionError(message)

    elif file.mode != 'rb':
        message = 'The file must be opened in binary mode instead of "{}".'.format(file.mode)
        raise TypeError(message)


    # ===================format check over==========================

    # b = base64.b64encode( s.encode('utf-8') ).decode('utf-8')
    if isinstance(file, str):
        f = open(file_str, 'rb')
        return hashlib.md5(f.read()).hexdigest()
    else:
        return hashlib.md5(file.read()).hexdigest()


# temp = open('/home/user/temp.txt', 'rb')
# r = md5(temp)
