import sys, os, time

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
CORE_DIR = os.path.join(THIS_DIR, '../../../')
sys.path.append(CORE_DIR)

from rsquarelabs_core.config import RSQ_DB_PATH, RSQ_SCRIPT_PATH, RSQ_BACKUP_PATH,\
    RSQ_PROJECTS_HOME, RSQ_EXPORT_PATH, RSQ_IMPORT_PATH, RSQ_HOME, RSQ_LOG_PATH
from rsquarelabs_core.engines.projects import Project
from rsquarelabs_core.engines.gromacs import Gromacs
from rsquarelabs_core.engines.db_engine import DBEngine
import unittest

# RSQ_DB_PATH = os.path.join(RSQ_HOME, 'rsquarelabs_test.db')
#
# db_object= DBEngine(RSQ_DB_PATH)

class BaseClass(unittest.TestCase):
    def setUp(self):
        print "-----> 1"
        self.project = Project()
        self.cur = DBEngine().cur

class ProjectTest(BaseClass):

    def test_project_create(self):
        print "-----> 2"

        super(ProjectTest, self).setUp()
        self.project.create(project_title='tak', project_tags='gpcr4', project_user_email='hello@world.com',
                            project_slug='tak', project_short_note='This is first project')

        print "testing create "

class RunTest(BaseClass):

    def test_run_create(self):
        print "-----> 4"
        super(RunTest, self).setUp()
        self.run_id = self.project.create_run(run_name="run_tak", version="1", parent_run_id="0",
                                master_id="0", run_data="This is data for running", is_delete="0",
                                project_id="1", protocol_class_name="ProteinMin")

        print "testing create of run"
        #return run_id


class GromacsTest(RunTest):

    def setUp(self):
        print "-----> 3"
        super(GromacsTest, self).setUp()
        super(GromacsTest, self).test_run_create()
        print self.run_id
        self.project_test = self.cur.execute("select id, path from projects").fetchone()
        print self.project_test
        #time.sleep(2)
        self.gromacs = Gromacs(project_id=self.project_test[0], working_dir=self.project_test[1], log_file=RSQ_LOG_PATH)

    def test_receptor_lig_min_mdp(self):
        print "-----> 5"
        self.gromacs.create_receptor_lig_min_mdp()
        print "creating receptor_lig_min_mdp"



def suite():
    test_suite = unittest.TestSuite()
    test_suite.addTest(unittest.makeSuite(ProjectTest))
    test_suite.addTest(unittest.makeSuite(GromacsTest))
    return test_suite

mySuit = suite()


runner = unittest.TextTestRunner()
runner.run(mySuit)
#
# if __name__ == '__main__':
#     unittest.main()
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
