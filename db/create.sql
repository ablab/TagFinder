create table t_status
    (scan_id int,
    protein_id int,
    evalue double precision,
    manual_score double precision,
    updated timestamp default now());

create table t_history
        (scan_id int,
        protein_id int,
        evalue double precision,
        version varchar(255),
        created timestamp default now());

create unique index idx_status on t_status(scan_id, protein_id);
create unique index idx_history on t_history(scan_id, protein_id, version);
