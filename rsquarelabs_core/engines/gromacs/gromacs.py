import shutil, argparse, sys, os
from rsquarelabs_core.utils import run_process
from core.messages import welcome_message, backup_folder_already_exists, \
    write_em_mpd_data, create_em_mdp_data
from core import settings
from rsquarelabs_core.engines.db_engine import DBEngine, RSQ_DB_PATH

"""
This module aimed at writing python wrapper around the tool Gromacs(www.gromacs.org) - a molecular dynamics package used
primarily in the field of Computational Drug Discovery also called as Computer Aided Drug Discovery(CADD) and others.

Brief about Gromacs :
GROMACS is a versatile package to perform molecular dynamics, i.e. simulate the Newtonian equations of motion for
systems with hundreds to millions of particles.It is primarily designed for biochemical molecules like proteins,
lipids and nucleic acids that have a lot of complicated bonded interactions.
Read more at http://www.gromacs.org/About_Gromacs


Read rsquarelabs-core/engines/gromacs/README.md for more info.
"""

db_object = DBEngine(RSQ_DB_PATH)

TOOL_NAME = 'r2_gromacs'

class ProteinLigMin(object):
    """
    ProteinLigandMin is a class that has all the methods to run a Protein-ligand Minimisation protocol using the tool -
    gromacs(gromacs.org).
    The module is designed based on tutorial by Justin Lemkul from Benvan Labs
    (http://www.bevanlab.biochem.vt.edu/Pages/Personal/justin/gmx-tutorials/complex/)
    """
    def __init__(self, *args, **kwargs):

        """
        To run the gromacs simulation of protein ligand complex, we need 3 files:
        1. protein structure file in .pdb format
        2. ligand structure file in *.gro format,
        3. ligand topology file in *.itp format
        :param args:
        :param kwargs:
        """

        self.ligand_file_path = kwargs.get('ligand_file', None)
        self.ligand_topology_file_path = kwargs.get('ligand_topology_file', None)
        self.protein_file_path = kwargs.get('protein_file', None)

        if self.ligand_file_path is None or self.ligand_topology_file_path is None or self.protein_file_path is None:
            print "Please provide 1. protein structure file, 2. ligand structure file, 3.ligand topology file"
            print "Terminating with input error "
            exit()

        """
        the working directory inside the project ie., if the project path is /home/me/rsquareProjects/moldyn and if 'run01'
        is given as working dir, it will create  /home/me/rsquareProjects/moldyn/run01 and will do all the runs, but by default
        all the jobs will executed in  /home/me/rsquareProjects/moldyn unless working_dir is specified
        """
        self.working_dir = kwargs.get('working_dir', './')
        self.verbose = kwargs.get('verbose', False)
        self.quiet = kwargs.get('quiet', False)
        self.project_id = kwargs.get('project_id', None)

        # A user cant use both the verbose and the quiet flag together
        if self.verbose is True and self.quiet is True:
            print 'Can\'t use both the verbose and quiet flags together'
            sys.exit()

    @staticmethod
    def welcome():
        """
        Prints out a welcome message, license info and the version.
        """
        print welcome_message

    @staticmethod
    def file_copy(source, destination):
        # TODO: There must be something better in the os module?
        in_file = open(source, 'r',0777)
        out_file = open(destination, 'w',0777)
        temp = in_file.read()
        out_file.write(temp)
        in_file.close()
        out_file.close()

    def gather_files(self):
        if not os.path.isfile(self.ligand_file_path):
            print 'Ligand file not found at ', self.ligand_file_path
            sys.exit()

        elif not os.path.isfile(self.ligand_topology_file_path):
            print 'Ligand Topology file not found at ', \
                self.ligand_topology_file_path
            sys.exit()

        elif not os.path.isfile(self.protein_file_path):
            print 'Protein file not found at ', self.protein_file_path
            sys.exit()

        else:
            print 'All data files found'

        if os.path.isdir(self.working_dir):
            print "Folder '" + self.working_dir + "' Aready exist"
            if os.path.isdir("BACKUP"):
                print "ERROR: Backup folder already exists :( "
                print backup_folder_already_exists
                sys.exit()
            else:
                if os.rename(self.working_dir, "BACKUP"):
                    print "Old " + self.working_dir + " was moved to BACKUP/"

        os.mkdir(self.working_dir)
        print "CHEERS: Working Directory " + self.working_dir + \
              " created Successfully"
        print "Moving the files to Working Directory" + self.working_dir
        shutil.copy2(self.protein_file_path, self.working_dir + 'protein.pdb')
        shutil.copy2(self.ligand_file_path, self.working_dir + 'ligand.gro')
        shutil.copy2(self.ligand_topology_file_path,
                     self.working_dir + 'ligand.itp')

    def create_topology(self):
        """
        This will create the topology for the protein, more about this command is explained in rsquarelabs-core/engines/gromacs/README.md

        Testcases:
        1. Takes protein.pdb as input and generates protein.gro, topol.top as output files

        Log file : step2.log

        :return:
        """
        print ">STEP1 : Initiating Procedure to generate topology for protein"
        pdb2gmx = settings.g_prefix + "pdb2gmx"
        step_no = "1"
        step_name = "Topology Generation"
        log_file =  self.working_dir + "step-%s.log"%step_no

        command = pdb2gmx + " -f " + self.working_dir + "protein.pdb -o " + \
            self.working_dir + "protein.gro -ignh -p " + \
            self.working_dir + "topol.top -i " + self.working_dir + \
            "posre.itp -ff gromos53a6 -water spc "
        print command
        print self.project_id
        run_process(step_no, step_name, command,TOOL_NAME, log_file, self.project_id)

    def prepare_system(self):
        """
        This will make the system that we will use for the rest of simulations

        1. This will combine both protein.gro + ligand.gro and makes system.gro
        2. includes ligand.itp into the topol,top
        3. Adds ligand identifier(typically UNK) into topol.top at the end

        Testcases:
        1. check if new file  'system.gro' is created
        2. Check for the presense of string '#include "ligand.itp"' in topol.top
        3. Check the ligand identifier(typically UNK ) presence in topol.top

        LogFile: None

        :return:
        """

        print ">STEP2 : Initiating Precedure to make system[Protein+Ligand]"
        start_from_line = 3  # or whatever line I need to jump to

        # TODO: WHAT IS THIS?
        protein = self.working_dir + "protein.gro"
        system = self.working_dir + "system.gro"
        ligand = self.working_dir + "ligand.gro"

        protein_file = open(protein, "r", 0777)
        ligand_file = open(ligand, "r",0777)
        system_file = open(system, 'wa', 0777)

        # get the last line of protein
        # get the count of Protein and Ligand files
        protien_lines_count = len(protein_file.readlines())
        ligand_lines_count = len(ligand_file.readlines())

        # print protien_lines_count
        # print ligand_lines_count
        # count of the system
        # TODO: Better name
        system_count = protien_lines_count + ligand_lines_count - 6
        protein_file.close()
        ligand_file.close()

        # open files for reading
        protein_file = open(protein, "r", 0777)
        ligand_file = open(ligand, "r", 0777)

        system_file.write(
            "System.gro Designed for Simulation by [bngromacs.py]\n")
        system_file.write(str(system_count) + "\n")

        line_counter = 1
        for line in protein_file:
            if line_counter in range(start_from_line,
                                     protien_lines_count):  # start_from_line :
                # print line
                system_file.write(line)
            line_counter += 1
        protein_file.close()

        line_counter = 1
        for line in ligand_file:
            if line_counter in range(start_from_line, ligand_lines_count):
                # print line
                system_file.write(line)
            line_counter += 1

            # get the last line of protein [the coordinates of the center]
        protein_file = open(protein, "r", 0777)
        last_line = protein_file.readlines()[-1]
        # print last_line
        system_file.write(last_line)
        print "CHEERS: system.gro WAS GENERATED SUCCESSFULLY"

        f1 = open(self.working_dir + 'topol.top', 'r',0777)
        f2 = open(self.working_dir + 'topol_temp.top', 'w',0777)
        for line in f1:
            f2.write(line.replace('; Include water topology',
                                  '; Include Ligand topology\n #include '
                                  '"ligand.itp"\n\n\n; Include water topology ')
                     )
        f1.close()
        f2.close()
        # swaping the files to get the original file
        f1 = open(self.working_dir + 'topol.top', 'w',0777)
        f2 = open(self.working_dir + 'topol_temp.top', 'r',0777)
        for line in f2:
            f1.write(line)
        f1.write("UNK        1\n")
        f1.close()
        f2.close()
        os.unlink(self.working_dir + 'topol_temp.top')
        print "INFO: Topology File Updated with Ligand topology info "
        print "CHEERS: STEP[2] SUCCESSFULLY COMPLETED :)\n\n\n"

    def solvate_complex(self):
        """
        This will define the unitcell for the simulation and solvates the cell/box.

        Testcases:
        1. creates files 'newbox.gro', 'solv.gro'

        Logfiles:
        1. step3.log
        2. step4.log

        :return:
        """


        print ">STEP3 : Initiating Procedure to Solvate Complex"
        editconf = settings.g_prefix + "editconf"
        step_no = "3"
        step_name = "Defining the Box"
        log_file = self.working_dir + "step-%s.log" % step_no
        command = editconf + " -f " + self.working_dir + "system.gro -o " + \
            self.working_dir + "newbox.gro -bt cubic -d 1 -c "
        run_process(step_no, step_name, command,TOOL_NAME, log_file, self.project_id)

        print ">STEP4 : Initiating Procedure to Solvate Complex"
        solvate = settings.g_prefix + "solvate"
        step_no = "4"
        step_name = "Solvating the Box"
        log_file = self.working_dir + "step-%s.log" % step_no
        command = solvate + " -cp " + self.working_dir + "newbox.gro -p " + \
            self.working_dir + "topol.top -cs spc216.gro -o " + \
            self.working_dir + "solv.gro   "

        run_process(step_no, step_name, command,TOOL_NAME, log_file, self.project_id)

    def write_em_mdp(self):
        """
        Writes the configuration file *.mdp needed for the run

        Testcase:
        1. file 'em.mdp' should be created

        :return:
        """

        print ">NOTE: Writing em.mdp"
        # TODO: Better name
        some_file = open(self.working_dir + "em.mdp", "w",0777)
        data = write_em_mpd_data
        some_file.write(str(data))
        some_file.close()



    def add_ions(self):
        """
        Add ions to the system to neutralise/add charge to the system.

        Testcase:
        1. creates output file 'ions.tpr'
        2. creates output file 'solv_ions.gro'
        
        :return:
        """

        print ">STEP5 : Initiating Procedure to Add Ions & Neutralise the " \
              "Complex"

        # TODO: Better name. Whats this?
        grompp = settings.g_prefix + "grompp"
        step_no = "5"
        step_name = "Check Ions "
        log_file = self.working_dir + "step-%s.log" % step_no

        # used to calc the charge, from the output
        temp_log_file = self.working_dir + "step-%s_temp.log" % step_no
        command = grompp + " -f " + self.working_dir + "em.mdp -c " + \
            self.working_dir + "solv.gro -p " + self.working_dir + \
            "topol.top -o " + self.working_dir + "ions.tpr -po " + \
            self.working_dir + "mdout.mdp > " + temp_log_file +" 2>&1"
        print command
        run_process(step_no, step_name, command,TOOL_NAME, log_file, self.project_id)

        # calculating the charge of the system
        # TODO: What is this doing? word??? Better name!
        word = 'total'  # Your word

        with open(temp_log_file,'r' ) as f:
            for line in f:
                if word in line:
                    s_line = line.strip().split()
                    two_words = (s_line[s_line.index(word) + 1],
                                 s_line[s_line.index(word) + 2])
                    charge = two_words[1]
                    break

        # TODO: This charge varibale might break the code
        print "Charge of the system is " + charge
        charge = float(charge)
        charge = int(round(charge))

        print "Charge =-=-=-=-=-=-=--=-=-=-==-%s" %charge

        if charge > 0:
            print "System has positive charge ."
            print "Adding " + str(charge) + " CL ions to Neutralize the system"
            genion = settings.g_prefix + "genion"
            step_no = "6"
            step_name = "Adding Negative Ions "
            log_file = self.working_dir + "step-%s.log" % step_no
            command = genion + " -s " + self.working_dir + "ions.tpr -o " + \
                self.working_dir + "solv_ions.gro -p " + self.working_dir + \
                "topol.top -nname CL -nn " + str(charge) +  " << EOF\nSOL\nEOF"
            run_process(step_no, step_name, command,TOOL_NAME, log_file, self.project_id)

        elif charge < 0:
            print "charge is negative"
            print "Adding " + str(-charge) + " CL ions to Neutralize the system"
            genion = settings.g_prefix + "genion"
            step_no = "6"
            step_name = "Adding Positive Ions "
            log_file = self.working_dir + "step-%s.log" % step_no
            command = genion + " -s " + self.working_dir + "ions.tpr -o " + \
                self.working_dir + "solv_ions.gro -p " + self.working_dir + \
                "topol.top -pname NA -np " + str(-charge) +  " << EOF\nSOL\nEOF"
            print command
            run_process(step_no, step_name, command,TOOL_NAME, log_file, self.project_id)

        elif charge == 0:
            print "System has Neutral charge , No adjustments Required :)"
            self.file_copy(self.working_dir +'ions.tpr', self.working_dir+ "solv_ions.tpr")

        print "DOUBLE CHEERS: SUCCESFULY PREPARED SYSTEM FOR SIMULATION"

    def write_emreal_mdp(self):
        # TODO: better name
        """
        Creates em_real.mdp which is used for minimisation

        Testcase:
        1. creates output file 'em_real.mdp'

        :return:
        """
        some_file = open(self.working_dir + "em_real.mdp", "w",0777)
        em_mdp = create_em_mdp_data
        some_file.write(em_mdp)
        print "CHEERS: em_real.mdp SUCCESSFULLY GENERATED :)"

    def minimize(self):
        """
        Minimises the protein-ligand system,

        Test cases:
        1. creates output file 'em.tpr'
        2. creates output files em.gro, em.edr , em.xtc, em.trr, em.log


        :return:
        """
        print "MESSAGE: "
        t = raw_input(
            "Did you check your complex !! do you wish to continue: (y/n)")

        if t == 'y':

            print ">STEP7 : Preparing the files for Minimisation"
            # grompp -f em_real.mdp -c solv_ions.gro -p topol.top -o em.tpr
            # mdrun -v -deffnm em
            grompp = settings.g_prefix + "grompp"
            mdrun = settings.g_prefix + "mdrun"
            step_no = "7"
            step_name = "Prepare files for Minimisation"
            log_file = self.working_dir + "step-%s.log" % step_no
            # max warn 3 only for now
            command = grompp + " -f " + self.working_dir + "em_real.mdp -c " +\
                self.working_dir + "solv_ions.gro -p " + self.working_dir\
                + "topol.top -o " + self.working_dir + "em.tpr -po " +\
                self.working_dir + "mdout.mdp -maxwarn 3 "
            print command
            run_process(step_no, step_name, command,TOOL_NAME, log_file, self.project_id)

            step_no = "8"
            step_name = " Minimisation"
            log_file = self.working_dir + "step-%s.log" % step_no

            command = mdrun + " -v  -s " + self.working_dir + "em.tpr -c " + \
                self.working_dir + "em.gro -o " + self.working_dir + \
                "em.trr -e " + self.working_dir + "em.edr -x " + \
                self.working_dir + "em.xtc -g " + self.working_dir + \
                "em.log "
            run_process(step_no, step_name, command,TOOL_NAME, log_file, self.project_id)
        else:
            print "Exiting on user request "
            sys.exit()


    """
    def nvt(self):
        print ">STEP9 : Initiating the Procedure to Equiliberate the System"
        print "Beginging Equiliberation with NVT Ensemble"
        grompp = settings.g_prefix + "grompp"
        mdrun = settings.g_prefix + "mdrun"
        step_no = "9"
        step_name = "Preparing files for NVT Equiliberation"
        command = grompp + "-f " + self.working_dir + "nvt.mdp -c " + \
            self.working_dir + "em.gro -p " + self.working_dir + \
            "topol.top -o " + self.working_dir + "nvt.tpr -po " + \
            self.working_dir + "mdout.mdp -maxwarn 3 > " + \
            self.working_dir + "step9.log 2>&1"
        run_process(step_no, step_name, command,TOOL_NAME, log_file, self.project_id)

        step_no = "10"
        step_name = "NVT Equiliberation"
        command = mdrun + " -v  -s " + self.working_dir + "nvt.tpr -c " + \
            self.working_dir + "nvt.gro -o " + self.working_dir + "nvt.trr -e "\
            + self.working_dir + "nvt.edr -x " + self.working_dir + \
            "nvt.xtc -g " + self.working_dir + "nvt.log > " + self.working_dir\
            + "step10.log 2>&1"
        run_process(step_no, step_name, command,TOOL_NAME, log_file, self.project_id)

    def npt(self):
        print ">STEP11 : Initiating the Procedure to Equiliberate the System"
        print "Beginging Equiliberation with NPT Ensemble"
        grompp = settings.g_prefix + "grompp"
        mdrun = settings.g_prefix + "mdrun"
        step_no = "11"
        step_name = "Preparing files for NPT Equiliberation"
        # grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr
        command = grompp + "-f " + self.working_dir + "npt.mdp -c " + \
            self.working_dir + "nvt.gro -p " + self.working_dir + \
            "topol.top -o " + self.working_dir + "npt.tpr -po " + \
            self.working_dir + "mdout.mdp -maxwarn 3 > " + self.working_dir + \
            "step11.log 2>&1"
        run_process(step_no, step_name, command,TOOL_NAME, log_file, self.project_id)

        step_no = "12"
        step_name = "NPT Equiliberation"
        command = mdrun + " -v  -s " + self.working_dir + "npt.tpr -c " + \
            self.working_dir + "npt.gro -o " + self.working_dir + \
            "npt.trr -e " + self.working_dir + "npt.edr -x " + \
            self.working_dir + "npt.xtc -g " + self.working_dir + "npt.log > "\
            + self.working_dir + "step12.log 2>&1"
        run_process(step_no, step_name, command,TOOL_NAME, log_file, self.project_id)

    def md(self):
        print "CHEERS :) WE ARE CLOSE TO SUCCESS :):)"
        print ">STEP13 : Initiating the Production Run"
        grompp = settings.g_prefix + "grompp"
        mdrun = settings.g_prefix + "mdrun"
        step_no = "13"
        step_name = "Preparing files for NPT Equiliberation"
        # grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr
        command = grompp + "-f " + self.working_dir + "md.mdp -c " + \
            self.working_dir + "npt.gro -p " + self.working_dir + \
            "topol.top -o " + self.working_dir + "md.tpr -po " + \
            self.working_dir + "mdout.mdp -maxwarn 3 > " + self.working_dir + \
            "step13.log 2>&1"
        run_process(step_no, step_name, command,TOOL_NAME, log_file, self.project_id)

        step_no = "14"
        step_name = "NPT Equiliberation"
        command = mdrun + " -v  -s " + self.working_dir + "md.tpr -c " + \
            self.working_dir + "md.gro -o " + self.working_dir + "md.trr -e " +\
            self.working_dir + "md.edr -x " + self.working_dir + "md.xtc -g " +\
            self.working_dir + "md.log > " + self.working_dir + \
            "step14.log 2>&1"
        run_process(step_no, step_name, command,TOOL_NAME, log_file, self.project_id)
    """

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-l', '--ligand',
                        help='Input a ligand file [*.gro]')
    parser.add_argument('-i', '--itp',
                        help='Input a ligand topology file [*.itp]')
    parser.add_argument('-p', '--protein',
                        help='Input a protein file (default:protein.pdb)')
    parser.add_argument('-w', '--wdir',
                        help='Working Directory of project (default:work)',
                        default='work/')
    parser.add_argument('-v', '--verbose', help='Loud and Noisy[default]',
                        action="store_true")
    parser.add_argument('-q', '--quiet', help='Be very quit',
                        action="store_true")

    arguments = parser.parse_args()

    """
      obj = ProteinLigMin(
        ligand_file=arguments.ligand,
        ligand_topology_file=arguments.itp,
        protein_file=arguments.protein,
        working_dir=arguments.wdir,
        verbose=arguments.verbose,
        quiet=arguments.quiet
    )

    obj.welcome()
    obj.gather_files()
    obj.create_topology()
    obj.prepare_system()
    obj.solvate_complex()
    obj.write_em_mdp()
    obj.add_ions()
    obj.write_emreal_mdp()
    obj.minimize()
    """


def hello():
    print "Hello World!!!"

def get_file_info(file):
    content = open(file).read()
    file_name = file.split("/")[-1]
    return [file_name, content]


def import_files(project_path, project_id):

    """
    This will import the files into the projects into the path. All the files imported via this method should be backed to
    somewhere in <HOME>/.rsquarelabs or the database.db, so that when the use clicks recreate project, that will be done.

    :param project_path: the project path to which the files should be imported
    :return:
    """

    protein_file_formats = [".pdb", ".gro"]
    ligand_file_formats = [".pdb", ".gro"]
    ligand_topology_file_formats = [".itp"]

    ligand_file = ""
    ligand_topology_file = ""
    protein_file = ""



    while not os.path.isfile(ligand_file):
        ligand_file = raw_input("Enter the path for ligand file : ")
        for format in ligand_file_formats:
            if ligand_file.endswith(format) and os.path.isfile(ligand_file):
                break

    while not os.path.isfile(ligand_topology_file):
        ligand_topology_file = raw_input("Enter the path for ligand topology file : ")
        for format in ligand_topology_file_formats:
            if ligand_file.endswith(format) and os.path.isfile(ligand_topology_file):
                break

    while not os.path.isfile(protein_file):
        protein_file = raw_input("Enter the path for protein file : ")
        for format in protein_file_formats:
            if protein_file.endswith(format) and os.path.isfile(protein_file):
                break

    # import this into the project
    ligand_structure_info = get_file_info(ligand_file)
    ligand_topology_info = get_file_info(ligand_topology_file)
    protien_structure_info = get_file_info(protein_file)


    # insert the file info ie., content and file_name

    db_object.do_insert("""
    INSERT INTO project_files(file_name, file_content, project_id)
    VALUES('%s','%s','%s')
    """ %(ligand_structure_info[0],ligand_structure_info[1],project_id))
    db_object.do_insert("""
    INSERT INTO project_files(file_name, file_content, project_id)
    VALUES('%s','%s','%s')
    """ %(ligand_topology_info[0],ligand_topology_info[1],project_id))

    db_object.do_insert("""
    INSERT INTO project_files(file_name, file_content, project_id)
    VALUES('%s','%s','%s')
    """ %(protien_structure_info[0],protien_structure_info[1],project_id))


    shutil.copy2(ligand_file, project_path)
    shutil.copy2(ligand_topology_file, project_path)
    shutil.copy2(protein_file, project_path)


    os.chmod(os.path.join(project_path, ligand_file), 0777)
    os.chmod(os.path.join(project_path, ligand_topology_file), 0777)
    os.chmod(os.path.join(project_path, protein_file), 0777)
