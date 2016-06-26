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

    def __init__(self, *args, **kwargs):
        print '----> 1'
        super(BaseClass, self).__init__(*args, **kwargs)
    #     #self.run_id = None
        self.project = Project(project_title='tak', project_tags='gpcr4', project_user_email='hello@world.com',
                            project_slug='tak', project_short_note='This is first project')
        self.cur = DBEngine().cur


    # def setUp(self):
    #     print "-----> 1"
    #     self.project = Project()
    #     self.cur = DBEngine().cur

class ProjectTest(BaseClass):

    def __init__(self, *args, **kwargs):
        super(ProjectTest, self).__init__(*args, **kwargs)

    # def setUp(self):
    #     self.run_id = None

    def test_project_create(self):
        print "-----> 2"
        self.project.save()

        print "testing create "

    def test_run_create(self):
        print "-----> 3"
        self.run_id = self.project.save_run(run_name="run_tak", version="1", parent_run_id="0",
                                master_id="0", run_data="This is data for running", is_delete="0",
                                project_id="1", protocol_class_name="ProteinMin")

        print "testing create of run"

    # def print_run_id(self):
    #     print self.run_id


class GromacsTest(BaseClass):

    # def __init__(self):

    def __init__(self, *args, **kwargs):
        print "-----> 4"
        super(GromacsTest, self).__init__(*args, **kwargs)
        # self.project_test = []
        self.run_id = None
        self.gromacs = None
        # super(GromacsTest, self).setUp()
        # self.projecttest_obj = ProjectTest()
        # self.run_id = self.projecttest_obj.run_id


    # def setUp(self):
    #     print '-----> set'
    #     print self._testMethodName
    #     self.project_test = self.cur.execute("select id, path from projects").fetchone()
    #     self.run_id = self.cur.execute("select run_id from runs").fetchone()[0]
    #     self.gromacs = Gromacs(project_id=self.project_test[0], receptor_file="/home/nitish/PycharmProjects/rsquarelabs-core/example_files/receptor.pdb",
    #                            working_dir=self.project_test[1], log_file=RSQ_LOG_PATH, run_id=self.run_id)

    # @classmethod
    # def setUpClass(cls):
    #     cls.select_project
    #     cls.select_run
    #     cls.gromacs_obj
    #
    # @classmethod
    def test_01_select_project(self):
        print "1"
        time.sleep(2)
        self.project_test = self.cur.execute("select id, path from projects").fetchone()
        print self.project_test

    def test_02_select_run(self):
        print "2"
        time.sleep(2)
        self.run_id = self.cur.execute("select run_id from runs").fetchone()[0]

    def test_03_gromacs_obj(self):
        print "3"
        time.sleep(2)

        print self.project_test
        self.gromacs = Gromacs(project_id=self.project_test[0],
                               receptor_file="/home/nitish/PycharmProjects/rsquarelabs-core/example_files/receptor.pdb",
                               working_dir=self.project_test[1], log_file=RSQ_LOG_PATH, run_id=self.run_id)

    # def test_project_create(self):
    #     print "-----> 2"
    #     self.project.create(project_title='tak', project_tags='gpcr4', project_user_email='hello@world.com',
    #                         project_slug='tak', project_short_note='This is first project')
    #
    #     print "testing create "
    #
    # def test_run_create(self):
    #     print "-----> 3"
    #     self.run_id = self.project.create_run(run_name="run_tak", version="1", parent_run_id="0",
    #                                           master_id="0", run_data="This is data for running", is_delete="0",
    #                                           project_id="1", protocol_class_name="ProteinMin")
    #
    #     print "testing create of run"

    # def test_something(self):
    #     print "-----> 6"
    #     self.project_id = "10"

    def test_04_receptor_lig_min_mdp(self):
        print "-----> 5"
        # print self.project_id
        self.gromacs.create_receptor_lig_min_mdp()
        print "creating receptor_lig_min_mdp"

    def test_05_gather_files(self):
        print "-----> 6"
        self.gromacs.gather_files()

    def test_06_pdb2gmax(self):
        print "-----> 7"
        input_name = "receptor.pdb"
        output_name = "receptor.gro"
        self.assertTrue(os.path.isfile(os.path.join(self.project_test[1], input_name)))
        self.gromacs.pdb2gmx(step_no=1, input_name=input_name, output_name=output_name, step_name="Creating topology for Protein",  parent_method_name="create_topology()", parent_method_serial=1)
        time.sleep(10)
        # self.assertTrue(os.path.isfile(os.path.join(self.project_test[1], output_name)))
        # os.path.isfile(self.)
        # self.assert
    #
    def test_04_editconf(self):
        time.sleep(10)
        print "-----> 8"
        input_name = "receptor.gro"
        output_name = "newbox.gro"
    #     self.assertTrue(os.path.isfile(os.path.join(self.project_test[1], input_name)))

        # self.gromacs.editconf(step_no=2, input_name=input_name, output_name=output_name, step_name="Defining the box",
        #               parent_method_name="create_water_box()", parent_method_serial=1)
        #
        # self.assertTrue(os.path.isfile(os.path.join(self.project_test[1], output_name)))

    # def test_05_solvate(self):
    #     time.sleep(10)
    #     print "-----> 9"
    #     input_name = "newbox.gro"
    #     output_name = "solv.gro"
    #     self.assertTrue(os.path.isfile(os.path.join(self.project_test[1], input_name)))
    #
    #     self.gromacs.solvate(step_no=3, input_name=input_name, output_name=output_name, step_name="Solvating the box",
    #                  parent_method_name="create_water_box()", parent_method_serial=2)
    #     self.assertTrue(os.path.isfile(os.path.join(self.project_test[1], output_name)))
    #
    # def test_06_grompp(self):
    #     time.sleep(10)
    #     print "-----> 10"
    #     input_name = "solv.gro"
    #     output_name = "ions.tpr"
    #     self.assertTrue(os.path.isfile(os.path.join(self.project_test[1], input_name)))
    #
    #     self.gromacs.grompp(step_no=4, input_name=input_name, output_name=output_name, mdp_file="ions.mdp",
    #                 step_name="Pre-processing to check the number of ions needed",
    #                 parent_method_name="neutralize_system()", parent_method_serial=1)
    #
    #     self.assertTrue(os.path.isfile(os.path.join(self.project_test[1], output_name)))
    #
    # def test_07_genion(self):
    #     time.sleep(10)
    #     print "-----> 11"
    #     input_name = "ions.tpr"
    #     output_name = "solv_ions.gro"
    #     self.assertTrue(os.path.isfile(os.path.join(self.project_test[1], input_name)))
    #
    #     self.gromacs.genion(step_no=5, input_name=input_name, output_name=output_name, step_name="Neutralizing the System",
    #                 parent_method_name="neutralize_system()", parent_method_serial=2)
    #     # self.grompp(step_no=6, input_name="solv_ions.gro", output_name="em.tpr", mdp_file="minim.mdp",
    #     #             step_name="Pre-processing the system before Minimisation", parent_method_name="minimize()",
    #     #             parent_method_serial=1)
    #
    #     self.assertTrue(os.path.isfile(os.path.join(self.project_test[1], output_name)))
    #
    # def test_08_mdrun(self):
    #     time.sleep(10)
    #     print "-----> 12"
    #     input_name = "em.tpr"
    #     self.assertTrue(os.path.isfile(os.path.join(self.project_test[1], input_name)))
    #
    #     self.gromacs.mdrun(step_no=7, input_name=input_name, nt=1, step_name="Final Minimisation",
    #                parent_method_name="minimize()", parent_method_serial=2)


def suite():
    suite = unittest.TestSuite()
    # suite.addTest(unittest.makeSuite(BaseClass))
    suite.addTest(unittest.makeSuite(ProjectTest))
    suite.addTest(unittest.makeSuite(GromacsTest))
    return suite

if __name__ == '__main__':
    runner = unittest.TextTestRunner(failfast=True)
    runner.run(suite())
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
