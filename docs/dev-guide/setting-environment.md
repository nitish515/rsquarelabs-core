# Development Env Setup


This is the first time setup doc for the developers.

## First time installation 

```bash
git clone https://github.com/rsquarelabs/core-client
virtualenv venv
source venv/bin/activate
# Use pip2.7 not pip, because pip might be linked to pip3 in some versions
pip2.7 install -r requirements/requirements.txt
pip2.7 install -r requirements/dev-requirements.txt
pip2.7 install -r requirements/testing-requirements.txt 

```

Install Gromacs (www.gromacs.org)

In Mac `brew install gromacs`

In ubuntu/debian `sudo apt-get install gromacs`

In Fedora/CentOS `sudo yum install gromacs`


## File system

- **rsquarelabs-core** - main module 
- **rsquarelabs-core/engines** - packages of the modules (gromacs, db_engine(sqlite)) 
- **rsquarelabs-core/utils** - utilities that can be used accross the core package
- **rsquarelabs-core/websuite** webserver that starts the webclient and can be accessed at http://localhost:9090
- **tests** - tests for all the modules/packages
- **sbin** - the python files in this module as supposed to be linked to /usr/local/bin/ so that they can be accessed from anywhere 



## Accessing the command

Once this package is installed via `pip install rsquarelabs-core` , all the python scripts in sbin are supposed to be accessbile from anywhere, 
but to simulate that behaviour in the development, we can create alias for these commands in `.bashrc` or add a softlink, here is how to
add these commands in `.bashrc` ie., `gedit ~/.bashrc`

add the following lines 

```bash
alias r2_gromacs='python2.7 /home/<USERNAME>/<R2-core-client-PATH>/sbin/r2_gromacs.py'
alias r2_server_start='python2.7 /home/<USERNAME>/<R2-core-client-PATH>/sbin/r2_server_start.py'
```

Example:

```bash
# In Mac
alias r2_gromacs='python2.7 /Users/rrmerugu/PycharmProjects/rsquarelabs-core/sbin/r2_gromacs.py'
alias r2_server_start='python2.7 /Users/rrmerugu/PycharmProjects/rsquarelabs-core/sbin/r2_server_start.py'

# In Linux
alias r2_gromacs='python2.7 /home/rrmerugu/PycharmProjects/rsquarelabs-core/sbin/r2_gromacs.py'
alias r2_server_start='python2.7 /home/rrmerugu/PycharmProjects/rsquarelabs-core/sbin/r2_server_start.py'

```

restart the terminal or do `source ~/.bashrc` to enable these commands for first time.


## Starting the client server 


Use the command `r2_server_start` to start the local client server, which is based on [bottle.py](http://bottlepy.org).
It is accessible at [http://localhost:9090](http://localhost:9090)




