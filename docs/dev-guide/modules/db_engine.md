#DBEngine Documentation

DBEngine is the memory(database) handler of this project i.e managing the whole databases through this module. In this
documentation we are going through the `DBEngine()` class. DBEngine deals with inserting, updating and selecting the data
in database.

#####Technologies used by application:
Python version - 2.7,
`sqlite3`

## Importing Modules

Importing some of the external modules like `sqlite3` and sub-module in `rsquarelabs_core.config` that are
needed in this module.

```python
import sqlite3, os, logging
from rsquarelabs_core.config import  RSQ_DB_LOG, RSQ_DB_PATH
```

## Provide Logging into this class

In order to track every execution in this class we need to provide logging before execution. Logging into log file is done
by using `debug()`, `info()` depending on execution.

```python
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
# create a file handler
handler = logging.FileHandler(RSQ_DB_LOG)
handler.setLevel(logging.INFO)
# create a logging format
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(handler)
```

## Creating object and Connection to the DataBase

Connecting to the database is done after by creating object to the `DBEngine()` class like `db_object = DBEngine(RSQ_DB_PATH)`.
Since connection to database needed a '.db' file i.e particular database file therefore `db_name` argument is needed to
this class. Here `RSQ_DB_PATH` is sent to the class as argument.

```python
def __init__(self, db_name=RSQ_DB_PATH):
    self.do_connect(db_name)
    """
    self.conn , self.cur are instalised in do_connect, but we are installising this here now,
    so we might need to change do_connect(), but all we are doing there is instalising
    """
    self.conn = sqlite3.connect(db_name)
    self.cur = self.conn.cursor()
```
This `self.cur` variable is created for manipulating the data which is called as `cursor()`.

Connection to the database is done by executing `db_object.do_connect(RSQ_DB_PATH)`. The table structure is design in
`schema.sql`.

```python
def do_connect(self, db_name):
    """
    This provides connection to the database.
    :param db_name: Path to the database file
    :return: Returns connection object to the database
    """
    self.conn = sqlite3.connect(db_name)
    sql_file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'schema.sql')
    tables_structure = open(sql_file_path).read()

    try:
        self.conn.executescript(tables_structure)
    except Exception as e:
        logger.debug(e)
        pass

    self.cur = self.conn.cursor()
    return self.conn
```

## Inserting Data into the Database

Inserting the given data into the database using the `cmd` and `param` parameters is done by executing `db_object.do_insert()`.
Returns the cursor.
```python
def do_insert(self, cmd, params):
    """
    This provides cursor for inserting the data into database using parameters.
    :param cmd: command (query) as a string for insertion.
    :param params: Constraint parameters for insertion data.
    :return: Returns cursor object which is inserted the data into database.
    """
    try:
        self.cur.execute(cmd, params)
        self.conn.commit()
    except Exception as e:
        logger.error(e)
        logger.error(e.message)
    return self.cur
```

**Note :** The `commit()` must be executed immediate after insertion in order to save the data into the database.

## Selecting Data from the DataBase

Selecting the data from the database using the `cmd` and `param` parameters is done by executing `db_object.do_select()`.
Returns the selected data.

```python
def do_select(self, cmd, params):
    """
    This provides data which is selected from database using given parameters.
    :param cmd: command (query) as a string for selecting.
    :param params: Constraint parameters for selecting data.
    :return: The selected data returned as a tuple.
    """
    data = self.cur.execute(cmd, params)
    return data
```
## Updating the Data into the DataBase

Updating the given data into the database using the `cmd` and `param` parameters is done by executing `db_object.do_update()`.
Returns the cursor.

```python
def do_update(self, cmd, params):
    """
    This provides cursor for updating the data into database using parameters.
    :param cmd: command (query) as a string for updating.
    :param params: Constraint parameters for updating data.
    :return: Returns cursor object which is updated the data into database.
    """
    try:
        self.cur.execute(cmd, params)
        self.conn.commit()
    except Exception as e:
        logger.error(e)
        logger.error(e.message)
    return self.cur
```

**Note :** The `commit()` must be executed immediate after updating in order to save the data into the database.