import sys, os

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
CORE_DIR = os.path.join(THIS_DIR, '../../../')
sys.path.append(CORE_DIR)

from rsquarelabs_core.config import RSQ_DB_PATH, RSQ_SCRIPT_PATH, RSQ_BACKUP_PATH,\
    RSQ_PROJECTS_HOME, RSQ_EXPORT_PATH, RSQ_IMPORT_PATH, RSQ_HOME, RSQ_LOG_PATH
from rsquarelabs_core.engines.projects import Project
from rsquarelabs_core.engines.gromacs import Gromacs
from rsquarelabs_core.engines.db_engine import DBEngine
import unittest, os

# RSQ_DB_PATH = os.path.join(RSQ_HOME, 'rsquarelabs_test.db')
#
# db_object= DBEngine(RSQ_DB_PATH)

class ProjectTest(unittest.TestCase):

    # def __init__(self):
    #     self.obj = ''
    #     self.cur = ''

    def setUp(self):
        self.obj = Project()
        self.cur = DBEngine().cur
        # self.proj = None

    def test_create(self):
        # self.objt = Project()
        # self.cur = self.obj.create(project_title='tak', project_tags='gpcr4', project_user_email='hello@world.com',
        #                     project_slug='tak', project_short_note='This is first project')

        print "testin create "


        # print "testing"
    # def test_select(self):
    #     self.proj = self.cur.execute("select id, path from projects").fetchone()

    # def tearDown(self):
    #     # print "+++++"
    #     os.remove(RSQ_DB_PATH)

# if __name__ == '__main__':
#     unittest.main()

class GromacsTest(ProjectTest, unittest.TestCase):

    def setUp(self):
        super(GromacsTest, self).setUp()
        # self.project = project
        # super(GromacsTest, self).test_create()
        self.project = self.cur.execute("select id, path from projects").fetchone()
        self.object = Gromacs(project_id=self.project[0], working_dir=self.project[1], log_file=RSQ_LOG_PATH)

    def test_gather_files(self):
        self.object.create_receptor_lig_min_mdp()

if __name__ == '__main__':
    unittest.main()
# unittest.main()
# #
#     def test_pdb2gmx():
#         // check if inout
#         obj.pdb2gmx()
#         // chec
# def something():
#     // check psd2gmx input if exist
#     // pdb2gmx
#     // check if pdb2gmx output if exist
