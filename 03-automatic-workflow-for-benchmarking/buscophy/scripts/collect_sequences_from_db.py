#!/usr/bin/env python3

import sqlite3
import os
import csv
import logging
import subprocess
import tempfile

def get_busco_duplication_stats(database):
    """
    Returns a dictionary: {Busco_id: duplicated_percentage}
    """
    with sqlite3.connect(database) as db_connection:
        cursor = db_connection.cursor()
        cursor.execute("""
            SELECT Busco_id,
                   COUNT(*) as total_count,
                   SUM(CASE WHEN Status = 'Duplicated' THEN 1 ELSE 0 END) as duplicated_count
            FROM full_table
            GROUP BY Busco_id
        """)
        stats = {}
        for busco_id, total_count, duplicated_count in cursor.fetchall():
            duplicated_percentage = (duplicated_count / total_count) * 100 if total_count else 0
            stats[busco_id] = duplicated_percentage
        return stats

def get_busco_ids_filtered(database, 
                                    min_taxa, 
                                    frag, 
                                    dupl, 
                                    min_genes,
                                    warn_callback)

    """
    This function filters species and BUSCO genes in three steps:

    1. Select species which fulfill the min_genes requirement 
       on (Complete or fragmented if frag or duplicated if dupl)
    2. Filter BUSCO IDs that have at least min_taxa among these species, 
    3. Filter BUSCO IDs further so that each species has at least min_genes BUSCO genes

    Returns:
    - final set of selected BUSCO IDs
    - set of valid species passing all criteria
    """
    with sqlite3.connect(database) as db_connection:
        cursor = db_connection.cursor()

        # 1. Identify species passing the min_genes threshold
        cursor.execute("SELECT DISTINCT Species FROM full_table")
        all_species = set(row[0] for row in cursor.fetchall())
        
        passing_species = set()
        for species in all_species:
            where_clause = ["Status='Complete'"]
            if frag:
                where_clause.append("Status='Fragmented'")
            if dupl:
                where_clause.append("Status='Duplicated'")

            query = (
                f"SELECT COUNT(DISTINCT Busco_id) FROM full_table WHERE Species=? "
                f"AND ({' OR '.join(where_clause)})"
            )
            cursor.execute(query, (species,))
            n_genes = cursor.fetchone()[0]
            
            if n_genes >= min_genes:
                passing_species.add(species)
            else:
                warn_callback(species, n_genes)
                logging.info(
                    f"WARNING: Species '{species}' has only {n_genes} qualifying BUSCO genes "
                    f"(min_genes={min_genes}). Excluded from output."
                )

        if not passing_species:
            return set(), set()


        
        species_placeholders = ','.join('?' for _ in passing_species)
        where_clause = ["Status='Complete'"]
        if frag:
            where_clause.append("Status='Fragmented'")
        if dupl:
            where_clause.append("Status='Duplicated'")

        query = (
            f"SELECT Busco_id, Species FROM full_table WHERE ({' OR '.join(where_clause)}) "
            f"AND Species IN ({species_placeholders})"
        )
        cursor.execute(query, tuple(passing_species))
        busco_to_species = {}
        for busco_id, species in cursor.fetchall():
            busco_to_species.setdefault(busco_id, set()).add(species)

        busco_ids = set()
        for busco_id, sps in busco_to_species.items():
            if len(sps) >= min_taxa:
                busco_ids.add(busco_id)

        if not busco_ids:
            return set(), set()

        return busco_ids, passing_species



def extract_from_db(database, selected_busco_ids, valid_species, frag, dupl):
    """
    Extract sequences from 'busco_sequences' table for the provided Busco_ids.
    Only sequences for valid_species are selected.
    Only sequences with Sequence_type 'amino' are selected.

    For duplicated BUSCOs, only the 'best_duplicate' (marked in the database column) is selected.
    DIAMOND is NOT invoked in this version.
    """
    if not selected_busco_ids or not valid_species:
        return []

    with sqlite3.connect(database) as db_connection:
        cursor = db_connection.cursor()

        placeholders_busco = ','.join(['?'] * len(selected_busco_ids))
        placeholders_species = ','.join(['?'] * len(valid_species))

        where_clause = ["Status='Complete'"]
        if frag:
            where_clause.append("Status='Fragmented'")
        if dupl:
            where_clause.append("Status='Duplicated'")

        # Get all relevant sequences, including best_duplicate info
        query = (
            f"SELECT b.Busco_id, b.Sequence, b.Species, b.best_duplicate "
            f"FROM busco_sequences b "
            f"INNER JOIN full_table f ON b.Busco_id = f.Busco_id AND b.Species = f.Species "
            f"AND ({' OR '.join(where_clause)}) "
            f"WHERE b.Sequence_type = 'amino' AND "
            f"b.Busco_id IN ({placeholders_busco}) AND "
            f"b.Species IN ({placeholders_species}) "
            f"ORDER BY b.Busco_id ASC;"
        )
        params = list(selected_busco_ids) + list(valid_species)
        cursor.execute(query, params)
        rows = cursor.fetchall()

    # Group rows by (Busco_id, Species)
    grouped = {}
    for busco_id, sequence, species, best_duplicate in rows:
        grouped.setdefault((busco_id, species), []).append((sequence, best_duplicate))

    final_results = []

    for (busco_id, species), seqs in grouped.items():
        if len(seqs) == 1:
            # No duplication - keep directly
            final_results.append((busco_id, seqs[0][0], species))
        else:
            # Choose only the best_duplicate marked as 1
            best_found = False
            for seq, best_duplicate in seqs:
                if best_duplicate == 1:
                    final_results.append((busco_id, seq, species))
                    best_found = True
                    break
            if not best_found:
                # Fallback: If no 'best_duplicate' is set, pick the first (should not occur if DB is marked correctly)
                final_results.append((busco_id, seqs[0][0], species))

    return final_results


def write_sequences_to_fasta(sequences, output_dir):
    """
    Writes sequences to separate FASTA files per Busco_id in the specified output directory.
    """
    os.makedirs(output_dir, exist_ok=True)
    current_busco = None
    fasta_out = None
    for busco_id, sequence, species in sequences:
        if busco_id != current_busco:
            if fasta_out:
                fasta_out.close()
            current_busco = busco_id
            destination = os.path.join(output_dir, f"{current_busco}.fas")
            fasta_out = open(destination, 'w')
        fasta_out.write(f">{species}\n{sequence}\n")
    if fasta_out:
        fasta_out.close()


def get_all_busco_ids(database, frag,dupl):
    """
    Returns a set of all Busco_ids in full_table that have at least one entry with
    Status 'Complete' (and 'Fragmented' if frag=True).
    """
    with sqlite3.connect(database) as db_connection:
        cursor = db_connection.cursor()
        query = (
            "SELECT DISTINCT Busco_id FROM full_table WHERE "
            "(Status == 'Complete' {fragments} {duplicates})"
        ).format(
            fragments="OR Status == 'Fragmented'" if frag else "",
            duplicates="OR Status == 'Duplicated'" if dupl else ""
        )
        cursor.execute(query)
        return set(row[0] for row in cursor.fetchall())

if __name__ == '__main__':
    database = snakemake.input.database
    frag = snakemake.params.fragmented
    dupl = snakemake.params.duplicated 
    min_taxa = int(snakemake.params.min_taxa)
    min_genes = int(snakemake.params.get("min_genes", 0))
    output_dir = snakemake.output[0]
    logfile = snakemake.log[0]


    logging.basicConfig(
        filename=logfile,
        level=logging.INFO,
        format="%(asctime)s %(levelname)s:%(message)s"
    )


    all_busco_ids = get_all_busco_ids(database, frag,dupl)

    warnings = []
    def warn_callback(species, n):
        msg = f"WARNING: Species '{species}' has only {n} qualifying BUSCO genes (min_genes={min_genes}). Excluded from output."
        print(msg)
        warnings.append(msg)

    selected_busco_ids, valid_species = get_busco_ids_filtered(
        database, 
        min_taxa, 
        frag, 
        dupl, 
        min_genes,
        warn_callback,
    )

    duplication_stats = get_busco_duplication_stats(database)

    with open(logfile, 'w') as f:
        print("Busco_id selection and duplication feedback:", file=f)
        print("------------------------------------------------", file=f)
        for busco_id in sorted(all_busco_ids):
            duplication = duplication_stats.get(busco_id, 0)
            warning = f"WARNING: {duplication:.2f}% Duplicated" if duplication > 30 else ""
            
            print(f"{busco_id}: SELECTED {warning}", file=f) if busco_id in selected_busco_ids else print(f"{busco_id}: NOT SELECTED {warning}", file=f)
        print("------------------------------------------------", file=f)
        print(f"Valid species: {sorted(valid_species)}", file=f)
        for w in warnings:
            print(w, file=f)


    sequences = extract_from_db(database, selected_busco_ids, valid_species , frag, dupl )
    
    write_sequences_to_fasta(sequences, output_dir)