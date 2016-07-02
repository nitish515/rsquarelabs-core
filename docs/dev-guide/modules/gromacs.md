<<<<<<< HEAD
# Gromacs Documentation

GROMACS is a software package for bio-molecules by performing very precise molecular dynamic simulations. We developed
every command related to gromacs into flexible methods in `Gromacs()` class in `init.py`. In this documentation we going
to discuss about Protein Minimisation `ProteinMin()` class in `gromacs.py`.

Every simulation steps has been recorded into log file.

#####Technologies used by application:
Python version - 2.7,
Gromacs version - 5

## Importing Modules

 Importing some of the external modules and sub-module in `rsquarelabs_core.utils`, `rsquarelabs_core.engines.db_engine`,
 `rsquarelabs_core.engines.gromacs` and `core.messages` that are needed in this module.

 ```python
import shutil, argparse, sys, os
from rsquarelabs_core.utils import run_process, get_file_info, import_files
from core.messages import welcome_message, backup_folder_already_exists, \
    write_em_mpd_data, create_em_mdp_data, ions_mdp, minim_mdp
from core import settings
from rsquarelabs_core.engines.db_engine import DBEngine, RSQ_DB_PATH
from rsquarelabs_core.engines.gromacs import Gromacs
import logging
 ```
Create the object for `DBEngine()` and defining the global variables.
```python
db_object = DBEngine(RSQ_DB_PATH)
TOOL_NAME = 'r2_gromacs'
```

## How to initiate Protein Minimisation

We initiate to `ProteinMin()` class creating object using needed arguments as of `Gromacs()` class because
`ProteinMin()` class has been inheritance from `Gromacs()`.

```python
obj = ProteinMin(project_id="PROJECT_ID", working_dir="WORKING_DIR", receptor_file="RECEPTOR_FILE",
                 log_file="LOG_FILE", run_id="RUN_ID")

```
####'__init__()' Operations in Gromacs class:

Checking the Working directory existence.

```python
if self.working_dir:
    if os.path.exists(self.working_dir):
        shutil.rmtree(self.working_dir)

    os.mkdir(self.working_dir,0777)
```
**Note :** Be careful using `shutil.rmtree()`.

Since every simulation steps has been recorded into log file. Here creating the object related to log file

```python
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
```

Gathering the files related to the run.

```python
def gather_files(self):
    if self.receptor_file_path:
        import_files(self.receptor_file_path, self.working_dir, self.project_id, self.run_id)
        os.chmod(self.working_dir, 0777)

```


## How to import protein file

We import the protein file i.e `receptor_file` into working directory by executing `obj.import_files()`.

```python
def import_files(self):
    """
    For this run we dont need many files, just protein is the input
    """
    protein_file_formats = ['.gro','.pdb']
    protein_file_path = self.receptor_file_path

    #Checking protein_file path exists as in given format or not. If does not exists, then enter the valid path.
    while not os.path.isfile(protein_file_path):
        protein_file_path = raw_input("Enter the path for protein file : ")
        for format in protein_file_formats:
            if protein_file_path.endswith(format) and os.path.isfile(protein_file_path):
                break
    import_files(protein_file_path, self.working_dir, self.project_id, self.run_id)
```

## How to Write ions.mdp and minim.mdp

We write ions.mdp and minim.mdp into working directory which are imported from `core.messages` executing
`obj.write_ions_mdp()` and `obj.write_prod_mdp()`.

```python
def write_ions_mdp(self):
    """
    Writes the configuration file *.mdp needed for the run

    Testcase:
    1. file 'ions.mdp' should be created

    """
    logging.info("Writing mdp files for protein minmisation ie., ions.mdp, minim.mdp")
    mdp_for_genion = open(os.path.join(self.working_dir, "ions.mdp"), "w", 0777)
    mdp_for_genion.write(str(ions_mdp))
    mdp_for_genion.close()

def write_prod_mdp(self):
    """
    Writes the configuration file *.mdp needed for the production run

    Testcase:
    1. file 'minim.mdp' should be created

    """
    mdp_for_min = open(os.path.join(self.working_dir, "minim.mdp"), "w", 0777)
    mdp_for_min.write(minim_mdp)
    mdp_for_min.close()
```
**Note :** We need to add run_and_record_process() to write_ions_mdp() and write_prod_mdp()

## Implementing Protein Minimisation Step-wise


### Creating Topology

Creating topology for the protein executing `obj.create_topology()` in order to use `pdb2gmx()` method inherited
from `Gromacs()` class.

```python
def create_topology(self):
        self.pdb2gmx(step_no=1, input_name="receptor.pdb", output_name="receptor.gro",
                     step_name="Creating topology for Protein",
                     parent_method_name="create_topology()", parent_method_serial=1)
```

The `pdb2gmx()` method inherited from `Gromacs()` class.
```python
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
        run_and_record_process( step_no, step_name, command, TOOL_NAME, log_file, self.project_id, self.run_id,
                                parent_method_name, parent_method_serial,command_method)

        self.logger.info("STEP%s: %s, completed. log written to %s " % (step_no, step_name, log_file))
```

### Creating Water Box

Creating water box for the providing water environment for further simulation. Its done by executing
`obj.create_water_box()` in order to use `editconf()` and `solvate()` methods inherited from `Gromacs()` class.

```python
def create_water_box(self):

        self.editconf(step_no=2, input_name="receptor.gro", output_name="newbox.gro",
                      step_name="Defining the box", parent_method_name="create_water_box()", parent_method_serial=1)
        self.solvate(step_no=3, input_name="newbox.gro", output_name="solv.gro",
                     step_name="Solvating the box", parent_method_name="create_water_box()", parent_method_serial=2)

```

The `editconf()` and `solvate()` methods inherited from `Gromacs()` class.

```python
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

    run_and_record_process( step_no, step_name, command, TOOL_NAME, log_file, self.project_id, self.run_id,
                            parent_method_name, parent_method_serial,command_method)

    self.logger.info("STEP%s: %s, completed. log written to %s " % (step_no, step_name, log_file))


def solvate(self, step_no, input_name="newbox.gro",
            output_name="solv.gro",
            step_name = "Solvating the Box",
            parent_method_name=None, parent_method_serial=1):

    set_file_premissions(os.path.join(self.working_dir,input_name))
    self.logger.info("STEP%s: Attempting the step %s " % (step_no, step_name))
    log_file = os.path.join(self.working_dir , "step-%s.log"%step_no)

    command_method = inspect.stack()[0][3]

    command = solvate + " -cp " + os.path.join(self.working_dir, input_name) + " -p " + \
              os.path.join(self.working_dir, "topol.top") + " -cs spc216.gro -o " + \
        os.path.join(self.working_dir , output_name)

    run_and_record_process( step_no, step_name, command, TOOL_NAME, log_file, self.project_id, self.run_id,
                           parent_method_name, parent_method_serial,command_method)

    self.logger.info("STEP%s: %s, completed. log written to %s " % (step_no, step_name, log_file))

```
###Creating Neutralize system

Creating the neutralize system by executing `obj.neutralize_system()`.

```python
def neutralize_system(self):
    self.write_ions_mdp()
    self.grompp(step_no=4, input_name="solv.gro", output_name="ions.tpr", mdp_file="ions.mdp",
                step_name="Pre-processing to check the number of ions needed", parent_method_name="neutralize_system()",
                parent_method_serial=1)
    self.genion(step_no=5, input_name="ions.tpr", output_name="solv_ions.gro", step_name="Neutralizing the System",
                parent_method_name="neutralize_system()", parent_method_serial=2)

```

The `grompp()` and `genion()` methods inherited from `Gromacs()` class.

```python

def genion(self, step_no, input_name="ions.tpr", output_name="solv_ions.gro",
            step_name = "Adding Ions to Neutralise the System", parent_method_name=None, parent_method_serial=1):

    set_file_premissions(os.path.join(self.working_dir,input_name))
    self.logger.info("STEP%s: Attempting the step %s " % (step_no, step_name))
    log_file = os.path.join(self.working_dir , "step-%s.log"%step_no)

    command_method = inspect.stack()[0][3]

    command = genion + " -s " + os.path.join(self.working_dir, input_name) + " -o " + \
              os.path.join(self.working_dir , output_name) + " -p " + \
              os.path.join(self.working_dir, "topol.top") + " -nname CL -pname NA -neutral << EOF\nSOL\nEOF"

    run_and_record_process( step_no, step_name, command, TOOL_NAME, log_file, self.project_id, self.run_id,
                            parent_method_name, parent_method_serial,command_method)

    self.logger.info("STEP%s: %s, completed. log written to %s " % (step_no, step_name, log_file))


def grompp(self,  step_no, input_name=None, output_name=None, mdp_file=None,
           step_name = "Gromacs Pre-processing", parent_method_name=None, parent_method_serial=1):

    set_file_premissions(os.path.join(self.working_dir,input_name))
    self.logger.info("STEP%s: Attempting the step %s " % (step_no, step_name))
    log_file = os.path.join(self.working_dir , "step-%s.log"%step_no)

    command_method = inspect.stack()[0][3]

    command = grompp + " -f " + os.path.join(self.working_dir, mdp_file) + " -c " + \
        os.path.join(self.working_dir, input_name) + " -p " + \
              os.path.join(self.working_dir , "topol.top")+ " -o " + os.path.join(self.working_dir , output_name) + \
              " -po " + os.path.join(self.working_dir, "mdout.mdp") + " -maxwarn 3"
    run_and_record_process( step_no, step_name, command, TOOL_NAME, log_file, self.project_id, self.run_id,
                            parent_method_name, parent_method_serial,command_method)

    self.logger.info("STEP%s: %s, completed. log written to %s " % (step_no, step_name, log_file))

```
### Minimisation

The final step minimisation is done by executing `obj.minimize()`.

```python
def minimize(self):
    self.write_prod_mdp()
    self.grompp(step_no=6, input_name="solv_ions.gro", output_name="em.tpr", mdp_file="minim.mdp",
                step_name="Pre-processing the system before Minimisation", parent_method_name="minimize()",
                parent_method_serial=1)
    self.mdrun(step_no=7, input_name="em.tpr", nt=1, step_name="Final Minimisation", parent_method_name="minimize()",
               parent_method_serial=2)

```
The `mdrun()` method inherited from `Gromacs()` class

```python
def mdrun(self, step_no, input_name=None, nt=1, step_name="mdrun ", parent_method_name=None, parent_method_serial=1):

    set_file_premissions(os.path.join(self.working_dir,input_name))
    input_with_no_extension = input_name.split(".")[0]

    self.logger.info("STEP%s: Attempting the step %s " % (step_no, step_name))
    log_file = os.path.join(self.working_dir , "step-%s.log"%step_no)

    command_method = inspect.stack()[0][3]

    command = mdrun + " -v  -s " + os.path.join(self.working_dir, input_with_no_extension) + ".tpr -c " + \
              os.path.join(self.working_dir , input_with_no_extension) +".gro -o " + \
              os.path.join(self.working_dir, input_with_no_extension) +".trr -e " + \
              os.path.join(self.working_dir , input_with_no_extension) +".edr -x " + \
              os.path.join(self.working_dir, input_with_no_extension) + ".xtc -g " + \
              os.path.join(self.working_dir, input_with_no_extension) + ".log  -nt " + str(nt)

    run_and_record_process( step_no, step_name, command, TOOL_NAME, log_file, self.project_id, self.run_id,
                            parent_method_name, parent_method_serial,command_method)

    self.logger.info("STEP%s: %s, completed. log written to %s " % (step_no, step_name, log_file))
```
=======
# Gromacs



## Protocol Sample 


For Gromacs 

```
import_files()
write_ions_mdp()
write_prod_mdp()
create_topology()
create_water_box()
neutralize_system()
minimize()
```
>>>>>>> 5a44461627866fb2e77352f45144cab3c3f82515
