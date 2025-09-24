#! /usr/bin/env Python3

##################
# IMPORT MODULES #
##################
import requests, tarfile, rich, os, zipfile, shutil, collections
import rich.progress
import rich_click as click
import subprocess as sp
from rich.progress import Progress, BarColumn, TextColumn

#############
# FUNCTIONS #
#############
def set_output_dir(output_):
    '''
    parses a user-provided string and returns the output directory
    '''
    try:
        if not output_.endswith('/'):
            output_ = output_ + '/'
        return f'{os.path.dirname(output_)}/'
    except AttributeError:
        return ''

def download_file(console, columns, url, output_directory, filename):
    '''
    Download a file from a given URL and save it to a local file.
    '''
    response = requests.get(url, stream=True)
    total_size = int(response.headers.get('content-length', 0))
    if len(url.split('/')[-1]) >= 5:
        terminal_filename = 'Downloading ' + url.split('/')[-1][0:5] + '...'
    else:
        spaces = ' ' * (8 - len(url.split('/')[-1]))
        terminal_filename = spaces + 'Downloading ' + url.split('/')[-1]
    try:
        with open(f'{output_directory}{filename}', 'wb') as file:
            with rich.progress.Progress(*columns) as progress_bar:
                task = progress_bar.add_task(console = console, description = f"[cyan]|{terminal_filename}[/] |", total=total_size)
                for chunk in response.iter_content(chunk_size=1024):
                    file.write(chunk)
                    progress_bar.update(task, advance=len(chunk))
    except FileNotFoundError as f:
        console.print(f"[cyan]|               ERROR[/] | [bold yellow]{f}, aborting analysis...[/]\n")
        exit()

def get_tar_file_count(tar_file):
    '''
    Count the number of files in the tarball
    '''
    with tarfile.open(tar_file, 'r:gz') as tar:
        return len([m for m in tar.getmembers() if m.isfile()])    

def tar_with_progress(console, columns, output_directory, tar_file):
    '''
    Extract tar file with progress bar
    '''
    total_files = get_tar_file_count(f'{output_directory}{tar_file}')
    if '/' in f'{output_directory}{tar_file}':
        command = ['tar', '-zxvf', f'{output_directory}{tar_file}', '-C', '/'.join(f'{output_directory}{tar_file}'.split('/')[:-1])]
    else:
        command = ['tar', '-zxvf', f'{output_directory}{tar_file}']
    with rich.progress.Progress(*columns) as progress_bar:
        task = progress_bar.add_task(console = console, description = "[cyan]|       Extracting...[/] |", total=total_files)
        with sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE, text=True) as proc:
            files_extracted = 0
            for line in proc.stderr:
                files_extracted += 1
                progress_bar.update(task, completed=files_extracted)

def remove_tar_intermediary(output_directory):
    '''
    remove intermediary files for the tar taxonomy file
    '''
    files_to_remove = ['citations.dmp', 'delnodes.dmp', 'division.dmp', 'gencode.dmp', 'merged.dmp', 'gc.prt', 'readme.txt', 'images.dmp']
    for file in files_to_remove:
        os.remove(f'{output_directory}{file}')

def blast_to_memory(console, columns, blast_input_):
    '''
    reads BLAST input file into memory and returns two dictionaries
    '''
    blast_dict = collections.defaultdict(dict)
    species_dict = collections.defaultdict(list)
    uniq_species_list = []
    total_size = os.path.getsize(blast_input_)
    with rich.progress.Progress(*columns) as progress_bar:
        pbar = progress_bar.add_task(console = console, description = "[cyan]|       Reading BLAST[/] |", total=total_size)
        with open(blast_input_, 'r') as infile:
            for line in infile:
                progress_bar.update(pbar, advance = len(line))
                if line == 'not assigned\n':
                    continue
                seq_id = line.split('\t')[0].split(';')[0]
                tax_id = line.split('\t')[2]
                pident = float(line.split('\t')[5])
                qcov = int(line.split('\t')[7])
                if seq_id not in blast_dict:
                    blast_dict[seq_id]['pident'] = pident
                    blast_dict[seq_id]['qcov'] = qcov
                    species_dict[seq_id].append(tax_id)
                    if tax_id not in uniq_species_list:
                        uniq_species_list.append(tax_id)
                elif pident == blast_dict[seq_id]['pident'] and qcov == blast_dict[seq_id]['qcov'] and tax_id not in species_dict[seq_id] and tax_id != 'N/A':
                    species_dict[seq_id].append(tax_id)
                    if tax_id not in uniq_species_list:
                        uniq_species_list.append(tax_id)
                elif pident > blast_dict[seq_id]['pident'] and qcov == blast_dict[seq_id]['qcov']:
                    blast_dict[seq_id]['pident'] = pident
                    blast_dict[seq_id]['qcov'] = qcov
                    species_dict[seq_id] = [tax_id]
                    if tax_id not in uniq_species_list:
                        uniq_species_list.append(tax_id)
                elif pident == blast_dict[seq_id]['pident'] and qcov > blast_dict[seq_id]['qcov']:
                    blast_dict[seq_id]['pident'] = pident
                    blast_dict[seq_id]['qcov'] = qcov
                    species_dict[seq_id] = [tax_id]
                    if tax_id not in uniq_species_list:
                        uniq_species_list.append(tax_id)
                elif pident > blast_dict[seq_id]['pident'] and qcov > blast_dict[seq_id]['qcov']:
                    blast_dict[seq_id]['pident'] = pident
                    blast_dict[seq_id]['qcov'] = qcov
                    species_dict[seq_id] = [tax_id]
                    if tax_id not in uniq_species_list:
                        uniq_species_list.append(tax_id)
    return blast_dict, species_dict, uniq_species_list

def fasta_to_memory(console, columns, fasta_):
    '''
    reads a fasta file of ASVs into memory and returns a list of sequence IDs
    '''
    zotu_list = []
    total_size = os.path.getsize(fasta_)
    with rich.progress.Progress(*columns) as progress_bar:
        pbar = progress_bar.add_task(console = console, description = "[cyan]|        Reading list[/] |", total=total_size)
        with open(fasta_, 'r') as infile:
            for line in infile:
                progress_bar.update(pbar, advance = len(line))
                if line.startswith('>'):
                    zotu_list.append(line.split(';')[0].lstrip('>'))
    return zotu_list

def table_to_memory(console, columns, table_input_):
    '''
    reads OTU/ASV table into memory and returns a list of sequence IDs
    '''
    zotu_list = []
    total_size = os.path.getsize(table_input_)
    with rich.progress.Progress(*columns) as progress_bar:
        pbar = progress_bar.add_task(console = console, description = "[cyan]|       Reading table[/] |", total=total_size)
        count = 0
        with open(table_input_, 'r') as infile:
            for line in infile:
                progress_bar.update(pbar, advance = len(line))
                count += 1
                if count == 1:
                    continue
                zotu_list.append(line.split('\t')[0])
    return zotu_list

def fill_out_blast_dict(console, columns, zotu_list, blast_dict, species_dict):
    '''
    adds entries to blast_dict based on OTU/ASV table list and returns updated blast_dict
    '''
    total_size = len(zotu_list)
    with rich.progress.Progress(*columns) as progress_bar:
        pbar = progress_bar.add_task(console = console, description = "[cyan]|   Fill missing data[/] |", total=total_size)
        for item in zotu_list:
            progress_bar.update(pbar, advance = 1)
            if item not in blast_dict:
                blast_dict[item]['pident'] = 'NA'
                blast_dict[item]['qcov'] = 'NA'
                species_dict[item].append('NA')
    return blast_dict

def crabs_to_memory(console, columns, crabs_, uniq_species_list):
    '''
    reads crabs_ to memory and extracts the necessary phylogenetic information
    '''
    species_lineages_dict = collections.defaultdict(dict)
    total_size = os.path.getsize(crabs_)
    with rich.progress.Progress(*columns) as progress_bar:
        pbar = progress_bar.add_task(console = console, description = "[cyan]|       Reading crabs[/] |", total=total_size)
        with open(crabs_, 'r') as infile:
            for line in infile:
                progress_bar.update(pbar, advance = len(line))
                tax_id = line.split('\t')[2]
                if tax_id in uniq_species_list:
                    if tax_id not in species_lineages_dict:
                        species_lineages_dict[tax_id]['kingdom'] = line.split('\t')[3]
                        species_lineages_dict[tax_id]['phylum'] = line.split('\t')[4]
                        species_lineages_dict[tax_id]['class'] = line.split('\t')[5]
                        species_lineages_dict[tax_id]['order'] = line.split('\t')[6]
                        species_lineages_dict[tax_id]['family'] = line.split('\t')[7]
                        species_lineages_dict[tax_id]['genus'] = line.split('\t')[8]
                        species_lineages_dict[tax_id]['species'] = line.split('\t')[9]
                        species_lineages_dict[tax_id]['common'] = line.split('\t')[1]
    return species_lineages_dict

def create_mrca(console, columns, species_dict, species_lineages_dict):
    '''
    looks up the taxonomic lineages for each sequence ID and returns the MRCA
    '''
    mrca_dict = {}
    contributing_species = collections.defaultdict(list)
    with rich.progress.Progress(*columns) as progress_bar:
        pbar = progress_bar.add_task(console = console, description = "[cyan]|    Calculating MRCA[/] |", total=len(species_dict))
        for item in species_dict:
            progress_bar.update(pbar, advance = 1)
            if len(species_dict[item]) == 1:
                lineage = species_lineages_dict[species_dict[item][0]]
                mrca_dict[item] = lineage
                try:
                    if lineage['species'] != 'NA':
                        contributing_species[item].append(lineage['species'])
                    else:
                        contributing_species[item].append(lineage['common'])
                except KeyError:
                    contributing_species[item].append('NA')
            else:
                count = 0
                mrca_lineage = []
                for tax_id in species_dict[item]:
                    lineage = species_lineages_dict[tax_id]
                    if lineage.get('species', 'NA') != 'NA':
                        contributing_species[item].append(lineage['species'])
                    else:
                        contributing_species[item].append(lineage.get('common', 'NA'))
                    count += 1
                    if count == 1:
                        mrca_lineage = lineage.copy()
                    else:
                        for key in mrca_lineage:
                            if key not in lineage or mrca_lineage[key] != lineage[key]:
                                mrca_lineage[key] = 'NA'
                mrca_dict[item] = mrca_lineage
    return mrca_dict, contributing_species

def write_output2(blast_dict, mrca_dict, contributing_species, output_):
    '''
    write results to output file
    '''
    separator = '\t'
    ranks_needed = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    with open(output_, 'w') as outfile:
        outfile.write(f'#OTU ID\t{separator.join(ranks_needed)}\tpident\tqcov\tmatching species IDs\n')
        for item in blast_dict:
            try:
                outfile.write(f'{item}\t{mrca_dict[item][ranks_needed[0]]}\t{mrca_dict[item][ranks_needed[1]]}\t{mrca_dict[item][ranks_needed[2]]}\t{mrca_dict[item][ranks_needed[3]]}\t{mrca_dict[item][ranks_needed[4]]}\t{mrca_dict[item][ranks_needed[5]]}\t{mrca_dict[item][ranks_needed[6]]}\t{blast_dict[item]["pident"]}\t{blast_dict[item]["qcov"]}\t{", ".join(contributing_species[item])}\n')
            except KeyError:
                outfile.write(f'{item}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n')

def names_to_memory(console, columns, names_input_):
    '''
    reads names.dmp to memory and returns two dictionaries
    '''
    id_key_dict = {}
    taxon_key_dict = {}
    total_size = os.path.getsize(names_input_)
    with rich.progress.Progress(*columns) as progress_bar:
        pbar = progress_bar.add_task(console = console, description = "[cyan]|       Extracting...[/] |", total=total_size)
        with open(names_input_, 'r') as infile:
            for line in infile:
                progress_bar.update(pbar, advance = len(line))
                if line.split('\t')[6] == 'scientific name':
                    tax_id = line.split('\t|\t')[0]
                    taxon_name = line.split('\t|\t')[1].replace(' ', '_')
                    id_key_dict[tax_id] = taxon_name
                    taxon_key_dict[taxon_name] = tax_id
    return id_key_dict, taxon_key_dict

def nodes_to_memory(console, columns, nodes_input_):
    '''
    reads nodes.dmp to memory and returns two dictionaries
    '''
    id_key_rank_up_values_dict = {}
    up_key_down_list_dict = collections.defaultdict(list)
    total_size = os.path.getsize(nodes_input_)
    with rich.progress.Progress(*columns) as progress_bar:
        pbar = progress_bar.add_task(console = console, description = "[cyan]|       Extracting...[/] |", total=total_size)
        with open(nodes_input_, 'r') as infile:
            for line in infile:
                progress_bar.update(pbar, advance = len(line))
                taxonId = line.split('\t|\t')[0]
                taxonIdUp = line.split('\t|\t')[1]
                taxonIdRank = line.split('\t|\t')[2]
                id_key_rank_up_values_dict[taxonId] = [taxonIdRank, taxonIdUp]
                up_key_down_list_dict[taxonIdUp].append(taxonId)
    return id_key_rank_up_values_dict, up_key_down_list_dict

def species_to_taxid_map(console, columns, uniq_species_list, taxon_key_dict):
    '''
    returns a mapping dict with species as key and tax ID as value
    '''
    species_taxid_dict = {}
    total_size = len(uniq_species_list)
    with rich.progress.Progress(*columns) as progress_bar:
        pbar = progress_bar.add_task(console = console, description = "[cyan]|       Extracting...[/] |", total=total_size)
        for species_name in uniq_species_list:
            progress_bar.update(pbar, advance = 1)
            if species_name.replace(' ', '_') in taxon_key_dict:
                species_taxid_dict[species_name] = taxon_key_dict[species_name.replace(' ', '_')]
            else:
                species_taxid_dict[species_name] = 'NA'
    return species_taxid_dict

def generate_lineage(console, columns, species_taxid_dict, id_key_rank_up_values_dict, id_key_dict):
    '''
    create the taxonomic lineage and return a dictionary
    '''
    species_lineages_dict = collections.defaultdict(list)
    total_size = len(species_taxid_dict)
    with rich.progress.Progress(*columns) as progress_bar:
        pbar = progress_bar.add_task(console = console, description = "[cyan]|       Extracting...[/] |", total=total_size)
        for item in species_taxid_dict.values():
            progress_bar.update(pbar, advance = 1)
            if item == 'NA':
                continue
            initial_item = item
            while item != id_key_rank_up_values_dict[item][1]:
                species_lineages_dict[initial_item].append([id_key_rank_up_values_dict[item][0], item, id_key_dict[item]])
                item = id_key_rank_up_values_dict[item][1]
    return species_lineages_dict

def filter_lineage(console, columns, species_lineages_dict):
    '''
    filters the existing taxonomic lineage and returns a clean dictionary
    '''
    filter_species_lineage_dict = collections.defaultdict(list)
    ranks_needed = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    total_size = len(species_lineages_dict)
    with rich.progress.Progress(*columns) as progress_bar:
        pbar = progress_bar.add_task(console = console, description = "[cyan]|       Extracting...[/] |", total=total_size)
        for lineage in species_lineages_dict:
            progress_bar.update(pbar, advance = 1)
            for rank in ranks_needed:
                rank_present = 0
                for item in species_lineages_dict[lineage]:
                    if rank == item[0]:
                        rank_present = 1
                        filter_species_lineage_dict[lineage].append(item[2])
                if rank_present == 0:
                    filter_species_lineage_dict[lineage].append('NA')
    return filter_species_lineage_dict

def generate_mrca(species_dict, filter_species_lineage_dict, species_taxid_dict):
    '''
    looks up the taxonomic lineages for each sequence ID and returns the MRCA
    '''
    mrca_dict = {}
    ranks_needed = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    for item in species_dict:
        if len(species_dict[item]) == 1:
            try:
                lineage = filter_species_lineage_dict[species_taxid_dict[species_dict[item][0]]]
            except KeyError:
                lineage = []
                for rank in ranks_needed:
                    lineage.append('NA')
            mrca_dict[item] = lineage
        else:
            count = 0
            mrca_lineage = []
            for species_name in species_dict[item]:
                lineage = filter_species_lineage_dict[species_taxid_dict[species_name]]
                if len(lineage) == 0:
                    continue
                count += 1
                if count == 1:
                    mrca_lineage = lineage
                else:
                    for i in range(len(mrca_lineage)):
                        if mrca_lineage[i] != lineage[i]:
                            mrca_lineage[i] = 'NA'
            mrca_dict[item] = mrca_lineage
    return mrca_dict

def write_output(blast_dict, mrca_dict, species_dict, output_):
    '''
    write results to output file
    '''
    separator = '\t'
    ranks_needed = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    with open(output_, 'w') as outfile:
        outfile.write(f'#OTU ID\t{separator.join(ranks_needed)}\tpident\tqcov\tmatching species IDs\n')
        for item in blast_dict:
            outfile.write(f'{item}\t{separator.join(mrca_dict[item])}\t{blast_dict[item]["pident"]}\t{blast_dict[item]["qcov"]}\t{", ".join(species_dict[item])}\n')