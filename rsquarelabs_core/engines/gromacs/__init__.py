import shutil, sys, os, inspect
from rsquarelabs_core.utils import run_process, run_and_record_process, get_file_info, import_files, \
    set_file_premissions
import logging
from core.messages import  write_em_mpd_data, create_em_mdp_data, ions_mdp, minim_mdp

# TODO-  Give option to change this via .yaml file in project
""" gmx for Gromacs 5.x , g_ for 4.5x """
g_prefix = "gmx "

TOOL_NAME = "r2_gromacs"

""" commands from the tool gromacs(www.gromacs.org) """
pdb2gmx     = g_prefix + "pdb2gmx"
editconf    = g_prefix + "editconf"
solvate     = g_prefix + "solvate"
grompp      = g_prefix + "grompp"
genion      = g_prefix + "genion"
mdrun       = g_prefix + "mdrun"



def check_inputs(input_files):
    logging.info("checking input_name")
    for file in input_files:
        if not os.path.exists(file):
            logging.error("File %s doesn't exist or no permissions to access the file " % file)
            print "File %s doesn't exist or no permissions to access the file " % file
            exit()


class Gromacs:
    """
    Gromacs parent class with gromacs commands as methods


    """
    def __init__(self, *args, **kwargs):
        """
        receptor and ligand might be changed into

        :param args:
        :param kwargs: receptor_file, ligand_file, ligand_topology_file, mdp_file_nvt,
                mdp_file_npt, mdp_file_min, mdp_file
        """

        self.receptor_file_path = kwargs.get('receptor_file', None)

        """ useful during receptor-ligand interaction study (https://en.wikipedia.org/wiki/Ligand_(biochemistry)) """
        self.ligand_file_path = kwargs.get('ligand_file', None)
        self.ligand_topology_file_path = kwargs.get('ligand_topology_file', None)

        """ mdp files required for respective steps """
        self.mdp_file_nvt = kwargs.get('mdp_file_nvt', None)
        self.mdp_file_npt = kwargs.get('mdp_file_npt', None)
        self.mdp_file_min = kwargs.get('mdp_file_min', None)
        self.mdp_file = kwargs.get('mdp_file', None)
        
        self.protocol_id = kwargs.get('protocol_id', None)

        """ Project directory / working directory """
        self.project_id = kwargs.get('project_id', None)
        self.working_dir = kwargs.get('working_dir', None)



        if self.working_dir:
            if os.path.exists(self.working_dir):
                shutil.rmtree(self.working_dir)
            print self.working_dir
            os.mkdir(self.working_dir,0777)

        # logger
        self.log_file = kwargs.get("log_file",None)

        if self.log_file:
            self.logger = logging.getLogger(__name__)
            self.logger.setLevel(logging.INFO)
            # create a file handler
            handler = logging.FileHandler(self.log_file)
            handler.setLevel(logging.INFO)
            # create a logging format
            formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            # add the handlers to the logger
            self.logger.addHandler(handler)

        # gather the files specified by the user
        self.gather_files()

        os.chmod(self.working_dir, 0777)







    def create_receptor_lig_min_mdp(self):

        self.logger.info("Writing mdp files for receptor-ligand ie., em.mdp, em_real.mdp")
        mdp_for_genion = open(self.working_dir + "em.mdp", "w", 0777)
        mdp_for_genion.write(str(write_em_mpd_data))
        mdp_for_genion.close()

        mdp_for_min = open(self.working_dir + "em_real.mdp", "w", 0777)
        mdp_for_min.write(create_em_mdp_data)
        mdp_for_min.close()

    #
    def gather_files(self):
        """
        This method is called during the init, if any files are parsed while initiaing the class, those files
        will be imported automatically.

        :return:
        """
        if self.receptor_file_path:
            import_files(self.receptor_file_path, self.working_dir, self.project_id, self.protocol_id)
            os.chmod(self.working_dir, 0777)

        if self.ligand_file_path:
            import_files(self.ligand_file_path, self.working_dir, self.project_id, self.protocol_id)

        if self.ligand_topology_file_path:
            import_files(self.ligand_topology_file_path, self.working_dir, self.project_id, self.protocol_id)



    def pdb2gmx(self, step_no,
                input_name="receptor.pdb",
                output_name="receptor.gro",
                step_name="Topology Generation", parent_method_name=None, parent_method_serial=1):

        set_file_premissions(os.path.join(self.working_dir,input_name))
        self.logger.info("STEP%s: Attempting the step %s " % (step_no, step_name))
        log_file =  os.path.join(self.working_dir , "step-%s.log"%step_no)

        command_method = inspect.stack()[0][3]

        command = pdb2gmx + " -f " + os.path.join(self.working_dir , input_name) + " -o " + \
            os.path.join(self.working_dir , output_name) + " -ignh -p " + \
                  os.path.join(self.working_dir, "topol.top" ) + " -i " +  \
                  os.path.join(self.working_dir ,"posre.itp") + " -ff gromos53a6 -water spc "
        run_and_record_process( step_no, step_name, command, TOOL_NAME, log_file, self.project_id, self.protocol_id, parent_method_name, parent_method_serial,command_method)

        self.logger.info("STEP%s: %s, completed. log written to %s " % (step_no, step_name, log_file))

    # def make_receptor_ligand_complex(self, step_no,
    #                                   receptor_name="receptor.gro",
    #                                   ligand_name="ligand.gro",
    #                                   step_name="Maing Receptor-Ligand Complex"):
    # 
    #     # set_file_premissions(os.path.join(self.working_dir,input_name))
    #     self.logger.info("STEP%s: Attempting the step %s " % (step_no, step_name))
    #     receptor = self.working_dir + receptor_name
    #     ligand = self.working_dir + ligand_name
    #     system = self.working_dir + "system.gro"
    # 
    # 
    #     protein_file = open(receptor, "r", 0777)
    #     ligand_file = open(ligand, "r",0777)
    #     system_file = open(system, 'wa', 0777)
    # 
    #     # get the last line of receptor
    #     # get the count of receptor and Ligand files
    #     protien_lines_count = len(protein_file.readlines())
    #     ligand_lines_count = len(ligand_file.readlines())
    # 
    # 
    #     # count of the system
    #     # TODO: Better name
    #     system_count = protien_lines_count + ligand_lines_count - 6
    #     protein_file.close()
    #     ligand_file.close()
    # 
    #     # open files for reading
    #     protein_file = open(receptor, "r", 0777)
    #     ligand_file = open(ligand, "r", 0777)
    # 
    #     system_file.write(
    #         "System.gro Designed for Simulation by [https://github.com/rsquarelabs/framework]\n")
    #     system_file.write(str(system_count) + "\n")
    # 
    #     start_from_line = 3  # or whatever line I need to jump to
    # 
    #     line_counter = 1
    #     for line in protein_file:
    #         if line_counter in range(start_from_line,
    #                                  protien_lines_count):  # start_from_line :
    #             # print line
    #             system_file.write(line)
    #         line_counter += 1
    #     protein_file.close()
    # 
    #     line_counter = 1
    #     for line in ligand_file:
    #         if line_counter in range(start_from_line, ligand_lines_count):
    #             # print line
    #             system_file.write(line)
    #         line_counter += 1
    # 
    #         # get the last line of receptor [the coordinates of the center]
    #     protein_file = open(receptor, "r", 0777)
    #     last_line = protein_file.readlines()[-1]
    #     # print last_line
    #     system_file.write(last_line)
    #     self.logger.info( "system.gro WAS GENERATED SUCCESSFULLY")
    # 
    #     f1 = open(self.working_dir + 'topol.top', 'r',0777)
    #     f2 = open(self.working_dir + 'topol_temp.top', 'w',0777)
    #     for line in f1:
    #         f2.write(line.replace('; Include water topology',
    #                               '; Include Ligand topology\n #include '
    #                               '"ligand.itp"\n\n\n; Include water topology ')
    #                  )
    #     f1.close()
    #     f2.close()
    #     self.logger.info("ligand topolgy added into topol.top ")
    #     # swaping the files to get the original file
    #     f1 = open(self.working_dir + 'topol.top', 'w',0777)
    #     f2 = open(self.working_dir + 'topol_temp.top', 'r',0777)
    #     for line in f2:
    #         f1.write(line)
    #     f1.write("UNK        1\n")
    #     f1.close()
    #     f2.close()
    #     os.unlink(self.working_dir + 'topol_temp.top')
    #     self.logger.info("STEP%s: %s, completed. log written to project log " % (step_no, step_name))
    # 

    def editconf(self, step_no,
                 input_name="system.gro",
                 output_name="newbox.gro",
                 step_name = "Defining the Box", parent_method_name=None, parent_method_serial=1):

        set_file_premissions(os.path.join(self.working_dir,input_name))
        
        self.logger.info("STEP%s: Attempting the step %s " %(step_no,step_name))
        log_file = os.path.join(self.working_dir , "step-%s.log"%step_no)

        command_method = inspect.stack()[0][3]

        command = editconf + " -f " + os.path.join(self.working_dir, input_name) + " -o " + \
            os.path.join(self.working_dir , output_name) + " -bt cubic -d 1 -c "

        print command
        run_and_record_process( step_no, step_name, command, TOOL_NAME, log_file, self.project_id, self.protocol_id, parent_method_name, parent_method_serial,command_method)
        self.logger.info("STEP%s: %s, completed. log written to %s " % (step_no, step_name, log_file))


    def solvate(self, step_no, input_name="newbox.gro", output_name="solv.gro", step_name = "Solvating the Box", parent_method_name=None, parent_method_serial=1):

        set_file_premissions(os.path.join(self.working_dir,input_name))
        self.logger.info("STEP%s: Attempting the step %s " % (step_no, step_name))
        log_file = os.path.join(self.working_dir , "step-%s.log"%step_no)

        command_method = inspect.stack()[0][3]

        command = solvate + " -cp " + os.path.join(self.working_dir, input_name) + " -p " + \
                  os.path.join(self.working_dir, "topol.top") + " -cs spc216.gro -o " + \
            os.path.join(self.working_dir , output_name)
        run_and_record_process( step_no, step_name, command, TOOL_NAME, log_file, self.project_id, self.protocol_id, parent_method_name, parent_method_serial,command_method)
        self.logger.info("STEP%s: %s, completed. log written to %s " % (step_no, step_name, log_file))


    def genion(self, step_no, input_name="ions.tpr", output_name="solv_ions.gro", step_name = "Adding Ions to Neutralise the System", parent_method_name=None, parent_method_serial=1):

        set_file_premissions(os.path.join(self.working_dir,input_name))
        self.logger.info("STEP%s: Attempting the step %s " % (step_no, step_name))
        log_file = os.path.join(self.working_dir , "step-%s.log"%step_no)

        command_method = inspect.stack()[0][3]

        command = genion + " -s " + os.path.join(self.working_dir, input_name) + " -o " + \
                  os.path.join(self.working_dir , output_name) + " -p " + \
                  os.path.join(self.working_dir, "topol.top") + " -nname CL -pname NA -neutral << EOF\nSOL\nEOF"
        print command
        run_and_record_process( step_no, step_name, command, TOOL_NAME, log_file, self.project_id, self.protocol_id, parent_method_name, parent_method_serial,command_method)
        self.logger.info("STEP%s: %s, completed. log written to %s " % (step_no, step_name, log_file))


    def grompp(self,  step_no, input_name=None, output_name=None, mdp_file=None,  step_name = "Gromacs Pre-processing", parent_method_name=None, parent_method_serial=1):

        set_file_premissions(os.path.join(self.working_dir,input_name))
        self.logger.info("STEP%s: Attempting the step %s " % (step_no, step_name))
        log_file = os.path.join(self.working_dir , "step-%s.log"%step_no)

        command_method = inspect.stack()[0][3]

        command = grompp + " -f " + os.path.join(self.working_dir, mdp_file) + " -c " + \
            os.path.join(self.working_dir, input_name) + " -p " + \
                  os.path.join(self.working_dir , "topol.top")+ " -o " + os.path.join(self.working_dir , output_name) + " -po " + \
                  os.path.join(self.working_dir, "mdout.mdp") + " -maxwarn 3"
        run_and_record_process( step_no, step_name, command, TOOL_NAME, log_file, self.project_id, self.protocol_id, parent_method_name, parent_method_serial,command_method)
        self.logger.info("STEP%s: %s, completed. log written to %s " % (step_no, step_name, log_file))


    def mdrun(self, step_no, input_name=None, nt=1, step_name="mdrun ", parent_method_name=None, parent_method_serial=1):

        set_file_premissions(os.path.join(self.working_dir,input_name))
        input_with_no_extension = input_name.split(".")[0]

        self.logger.info("STEP%s: Attempting the step %s " % (step_no, step_name))
        log_file = os.path.join(self.working_dir , "step-%s.log"%step_no)

        command_method = inspect.stack()[0][3]

        command = mdrun + " -v  -s " + os.path.join(self.working_dir, input_with_no_extension) + ".tpr -c " + \
                  os.path.join(self.working_dir , input_with_no_extension) +".gro -o " + \
                  os.path.join(self.working_dir, input_with_no_extension) +".trr -e " + os.path.join(self.working_dir , input_with_no_extension) +".edr -x " + \
                  os.path.join(self.working_dir, input_with_no_extension) + ".xtc -g " + \
                  os.path.join(self.working_dir, input_with_no_extension) + ".log  -nt " + str(nt)
        run_and_record_process( step_no, step_name, command, TOOL_NAME, log_file, self.project_id, self.protocol_id, parent_method_name, parent_method_serial,command_method)
        self.logger.info("STEP%s: %s, completed. log written to %s " % (step_no, step_name, log_file))