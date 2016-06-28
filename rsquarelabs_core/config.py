__author__ = 'rrmerugu'
import os

USER_HOME_FOLDER = os.getenv('HOME')

RSQ_PROJECTS_HOME = os.path.join(USER_HOME_FOLDER, 'rsquarelabsProjects')
RSQ_PROJECTS_CONFIG = os.path.join(RSQ_PROJECTS_HOME, '.config.json') # not very much needed
RSQ_HOME = os.path.join(USER_HOME_FOLDER, '.rsquarelabs')
RSQ_LOG_PATH = os.path.join(RSQ_HOME, 'rsquarelabs.log')

RSQ_BACKUP_PATH = os.path.join(RSQ_HOME, 'backup')
RSQ_SCRIPT_PATH = os.path.join(RSQ_HOME, 'automator_scripts')
RSQ_DB_PATH = os.path.join(RSQ_HOME, 'rsquarelabs.db')
RSQ_DB_LOG = os.path.join(RSQ_HOME, 'DBEngine.log')

RSQ_EXPORT_PATH = os.path.join(RSQ_HOME, 'export')
RSQ_IMPORT_PATH = os.path.join(RSQ_HOME, 'import')


HUB_API_BASE_URL = 'http://localhost:8000/'
HUB_API_BASE_URL.rstrip('/')