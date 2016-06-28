CREATE TABLE projects
       (id INTEGER PRIMARY KEY     AUTOINCREMENT,
       title          TEXT    NOT NULL,
       tags           TEXT    NOT NULL,
       user_email     TEXT    NOT NULL,
       short_note     TEXT    NOT NULL,
       slug           TEXT    NOT NULL,
       path           TEXT,
       config         TEXT,
       log            TEXT,
       type           TEXT,
       date           TEXT,
       is_delete      INT);
CREATE TABLE project_activity
    ( id  INTEGER PRIMARY KEY     AUTOINCREMENT,
    tool_name TEXT NOT NULL,
    step_no TEXT NOT NULL,
    step_name TEXT NOT NULL,
    command TEXT NOT NULL,
    pid TEXT,
    status TEXT NOT NULL,
    log_file TEXT NOT NULL,
    project_id INTEGER NOT NULL,
    created_at TEXT,
    updated_at TEXT,
    pid_status TEXT,
    run_id INT NOT NULL,
    parent_method_name TEXT NOT NULL,
    parent_method_serial INT NOT NULL,
    command_method TEXT NOT NULL);
CREATE TABLE project_files
    (id INTEGER PRIMARY KEY AUTOINCREMENT,
    file_name TEXT NOT NULL,
    file_content BLOB NOT NULL,
    project_id INTEGER NOT NULL,
    created_at TEXT,
    updated_at TEXT,
    run_id INT NOT NULL);
CREATE TABLE runs
    (run_id INTEGER PRIMARY KEY AUTOINCREMENT,
    run_name TEXT NOT NULL,
    version INT,
    parent_run_id INT,
    master_id INT,
    run_data TEXT NOT NULL,
    is_delete INT,
    python_file TEXT,
    log_file TEXT,
    w_dir TEXT,
    project_id INT NOT NULL,
    class_name TEXT NOT NULL,
    run_description TEXT );
CREATE TABLE protocols
    (protocol_id INTEGER PRIMARY KEY AUTOINCREMENT,
    protocol_name TEXT,
    class TEXT,
    protocol_data TEXT,
    description TEXT,
    author TEXT,
    date TEXT);
CREATE TABLE notes
    (note_id INTEGER PRIMARY KEY AUTOINCREMENT,
    note TEXT,
    run_id INT NOT NULL)
