# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 17:50:28 2018

@author: pedro
"""

#%%  
from Bio import Entrez #Para aceder ao NCBI
from Bio import SeqIO
from Bio import SeqFeature
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW
#%%
#pesquisa automatizada de artigos no Pubmed de interesse
Entrez.email = "...@aaa.org"
handle = Entrez.efetch(db="pubmed", id="28588310,29111265,23519164", retmode="xml")
search_n315=Entrez.read(handle)

#%%
#genes essenciais presentes na simulação de genes essenciais do modleo do repositórii do Optflux
# com função objetivo R_biomass_SA_8a 
genes_essenciais=['SA0923','SA0924','SA0925','SA0926','SA0920','SA0921','SA0922','SA0938','SA0937','SA0916',
                  'SA0917','SA0918','SA0919','SA0912','SA0913','SA0915','SA0910','SA0911','SA1938','SA0842',
                  'SA0843','SA0965','SA1297','SA1299','SA1298','SA1177','SA1052','SA2027','SA1065','SA1397',
                  'SA1150','SA2127','SA1165','SA1164','SA1166','SA2136','SA2186','SA2287','SA1197','SA1199',
                  'SA2288','SA1088','SA1089','SA1571','SA0486','SA0244','SA0487','SA1461','SA1585','SA0375',
                  'SA1228','SA1348','SA1227','SA1229','SA1587','SA0376','SA0134','SA1586','SA1347','SA1226',
                  'SA1589','SA0016','SA1588','SA1104','SA1346','SA1558','SA0347','SA1439','SA2406','SA0344',
                  'SA0345','SA0346','SA0592','SA0472','SA0593','SA0473','SA0594','SA0474','SA1205','SA1202',
                  'SA2412','SA0596','SA1201','SA0597','SA1204','SA1203','SA2413','SA1494','SA2341','SA2465',
                  'SA1496','SA2464','SA1493','SA1250','SA1492','SA2467','SA2466','SA1259','SA2468','SA2347',
                  'SA0176','SA0177','SA2470','SA2471','SA0178','SA0179','SA1352','SA1115','SA1244','SA2333',
                  'SA2456','SA1245','SA1487','SA2334','SA1126','SA1368','SA1858','SA1735','SA1731','SA1860',
                  'SA1749','SA0419','SA1865','SA1861','SA1864','SA1863','SA0506','SA1959','SA0865','SA1608',
                  'SA1729','SA1728','SA0512','SA1965','SA0996','SA0997','SA0756','SA1724','SA0994','SA0995',
                  'SA1651','SA0683','SA1650','SA1652','SA1412','SA0693','SA1309','SA0457','SA0458','SA1427',
                  'SA1669','SA1306','SA1301','SA1422','SA1545','SA1424','SA0549','SA0547','SA0669','SA0548',
                  'SA0670','SA0793','SA0439','SA0794','SA0795','SA1523','SA0796','SA1522']#optflux
genes_essenciais=sorted(genes_essenciais) #ordena alfabeticamente os genes
Entrez.email='lol@lul.com'
handle=Entrez.efetch(db='nucleotide', rettype='gb', retmode='text', id='BA000018.3') #retira o genoma do S.aureus N315
aureus=SeqIO.read(handle, 'genbank') #lê o genoma
SeqIO.write(aureus,"saureusn315.gbk","genbank") #escreve num ficheiro o genoma
handle.close()
#%%
#imprime dados de anotação do genoma
print(aureus.id)
print(aureus.name)
print()
print(aureus.description)
print()
print(aureus.annotations)
print()
print(aureus.dbxrefs)

#%%
#guarda as regiões codificantes dos genes essenciais do optflux, através da informações dos qualifiers dos features de cada CDS 
pos_genes_essenciais=[]
genesaureus=[]
for i in range(len(aureus.features)):
    if aureus.features[i].type=='CDS':
        essencial=False
        pos=0
        for g_essencial in genes_essenciais:
            if g_essencial in str(aureus.features[i].qualifiers["note"]):
                essencial=True
                pos=genes_essenciais.index(g_essencial)
                break
        if essencial:
            genesaureus.append(i)
            pos_genes_essenciais.append(pos)
genes=[]
for i in genesaureus:
    genes.append(aureus.features[i].extract(aureus.seq))

#%%
#para cada CDS imprime os qualifires
for i in pos_genes_essenciais:
    print(aureus.features[i])
    print(aureus.features[i].qualifiers)

#%%
E_VALUE_THRESH = 0.05 # assume-se que o threshold sobre o qual a homologia é insignificativa é 0.05
blast_records_BA=[]
#%%
#para cada gene essencial do Optflux, faz-se o Blast contra o genoma humano (txid9606 [ORGN]) e guarda o alinhamento
for gene in genes:
    result_handle= NCBIWWW.qblast("blastn","nr",gene._data,entrez_query='txid9606 [ORGN]', hitlist_size=25)
    blast_records=NCBIXML.read(result_handle)
    try:
        print ('Homologia.')
        blast_records_BA.append((blast_records.alignments))
    except:
        print('Blast não possível.')
        blast_records_BA.append('Blast não possível.')
#%%      
# verifica-se se os alinhamentos são significativos (com valor de evalue/expect menor que 0.05)
genes_sem_homologia=[]
for i in range(len(blast_records_BA)):
    print(genes_essenciais[i], end= " ")
    try:
       print(blast_records_BA[i][0].hsps[0].expect)
       if blast_records_BA[i][0].hsps[0].expect>E_VALUE_THRESH:
           genes_sem_homologia.append([genes_essenciais[i],blast_records_BA[i][0].hsps[0].expect]) #i+84
    except:
        print("SEM RESULTADOS.")

#%%
#verifica-se se os genes essenciais obtidos com o Optflux batem certo com os presentes na BD DEG 
genes_DEG=open("genes_essenciais_DEG.txt","r")
genes_essenciais_deg=genes_DEG.readlines()
genes_DEG.close()
g_essenciais_deg=[]
for i in range(len(genes_essenciais_deg)):
    if genes_essenciais_deg[i].find("/")==-1:
        genes_essenciais_deg[i]=genes_essenciais_deg[i].strip('\n')
        g_essenciais_deg.append(genes_essenciais_deg[i])
    else:
        genes_essenciais_deg[i]=genes_essenciais_deg[i].strip('\n')
        temp=genes_essenciais_deg[i].split("/")
        for j in temp:
            g_essenciais_deg.append(j)
genes_essenciais_deg=g_essenciais_deg

#%%
genes_optflux_gb=[]
for i in range(len(genes_essenciais_deg)):
    for j in range(len(aureus.features)):
        if aureus.features[j].type=="CDS":
            if (str(genes_essenciais_deg[i]) in str(aureus.features[j].qualifiers["note"]) or
                str(genes_essenciais_deg[i]) in str(aureus.features[j].qualifiers["gene"])):
                    genes_optflux_gb.append([genes_essenciais_deg[i],j])
                    break
#%%
genes_deg=[]
for i in genes_optflux_gb:
    ind = i[1]
    string=(str(aureus.features[ind].qualifiers["note"]))
    ind_2=string.find("ORFID:")
    genes_deg.append(string[ind_2+6:ind_2+12])

#%%
deg_optflux=[]
genes_s_homologia_optflux=[]
#for i in result_s_homologia: genes_s_homologia_optflux.append(i[0])
for i in genes_sem_homologia: genes_s_homologia_optflux.append(i[0])
genes_s_homologia_optflux=sorted(genes_s_homologia_optflux)
for i in genes_deg:
    if i in genes_s_homologia_optflux:deg_optflux.append(i)
    
#%%
#   Genes sem homologia com humanos, presentes no Optflux e DGE
# 'SA0179', 0.309651  https://www.uniprot.org/uniprot/P60296
# 'SA0457', 0.267056 https://www.uniprot.org/uniprot/Q7A7B4
# 'SA0924', 0.14368 https://www.uniprot.org/uniprot/P99162
# 'SA0997', 0.638148 https://www.uniprot.org/uniprot/P63638
# 'SA1104', 0.205002 https://www.uniprot.org/uniprot/Q7A5Y4
# 'SA1177', 1.87107 https://www.uniprot.org/uniprot/P99161
    
#['SA1204', 0.281624 https://www.uniprot.org/uniprot/P66987
# 'SA1259', 0.372101 https://www.uniprot.org/uniprot/P99079
# 'SA1492', 0.889915 https://www.uniprot.org/uniprot/P64334
# 'SA1494', 0.244495 https://www.uniprot.org/uniprot/P64341
# 'SA1728', 0.429948 https://www.uniprot.org/uniprot/P99150
# 'SA2027', 0.588693 https://www.uniprot.org/uniprot/P99062
# 'SA2406'] 0.391652 https://www.uniprot.org/uniprot/A0A0H3JNT1
    #########
# 'SA1522', 0.0679558    
# 'SA1346', 0.0967526

#%%
# lê o ficheiro do Proteome do S. aureus N315
aureus_prot=SeqIO.parse("uniprot-proteome_UP000000751.fasta", "fasta")
#%%
aa_fasta=[]
for i in aureus_prot:
    aa_fasta.append(i)


#%%
#retira os dados das proteínas do Uniprot
#https://www.uniprot.org/uniprot/P99079
#SA1259  Dihydrofolate reductase  
# annotation score: 3/5

 
from Bio import ExPASy
from Bio import SwissProt
handle = ExPASy.get_sprot_raw('P99079')
record = SwissProt.read(handle)

#print(record.description)
#for ref in record.references:
#    print("authors:", ref.authors)
#    print("title:", ref.title)
#print(record.organism_classification)

#dir(record)

#%%
info=[record.accessions, record.annotation_update, record.comments, record.created, record.cross_references, record.data_class, record.description, record.entry_name, record.features, record.gene_name, record.host_organism, record.host_taxonomy_id, record.keywords, record.molecule_type, record.organelle, record.organism, record.organism_classification, record.protein_existence, record.references, record.seqinfo, record.sequence, record.sequence_length, record.sequence_update, record.taxonomy_id]
for i in range(len(info)):
    print(dir(record)[i+26],":" ,end=" ")
    if info[i] != record.references: 
        print(info[i])
        print()
    else:
        for k in range(len(info[i])):
            print()
            print("Authors:" ,end=" ")
            print(info[i][k].authors)
            print("comments:" ,end=" ")
            print(info[i][k].comments)
            print("location:" ,end=" ")
            print(info[i][k].location)
            print("number: ", end=" ")
            print(info[i][k].number)
            print("positions:" ,end=" ")
            print(info[i][k].positions)
            print("references:" ,end=" ")
            print(info[i][k].references)
            print("title:" ,end=" ")
            print(info[i][k].title)
            print()

#%%
#retira os dados do NCBI CDD para essa proteína
handle=Entrez.efetch(db='protein', rettype='gb', retmode='text', id='P99079')
Dhf_reductase=SeqIO.read(handle, 'genbank')
SeqIO.write(Dhf_reductase,"Dihydrofolate_reductase.gbk","genbank")
handle.close()
for i in Dhf_reductase.features:
    print(i)
    print(i.qualifiers)
    
    
#%%    
#prot_ncbi= [Dhf_reductase.annotations,
# Dhf_reductase.dbxrefs,
# Dhf_reductase.description,
# Dhf_reductase.features,
# Dhf_reductase.format,
# Dhf_reductase.id,
# Dhf_reductase.letter_annotations,
# Dhf_reductase.lower,
# Dhf_reductase.name,
# Dhf_reductase.reverse_complement,
# Dhf_reductase.seq,
# Dhf_reductase.translate,
# Dhf_reductase.upper]   

    
#%%
#retira os dados das proteínas do Uniprot
# SA0997 Glutamate racemase
# https://www.uniprot.org/uniprot/P63638
# annotation score: 2/5
handle = ExPASy.get_sprot_raw('P63638')
record = SwissProt.read(handle)
#print(record.description)
#for ref in record.references:
#    print("authors:", ref.authors)
#    print("title:", ref.title)
#print(record.organism_classification)

#%%
info=[record.accessions, record.annotation_update, record.comments, record.created, record.cross_references, record.data_class, record.description, record.entry_name, record.features, record.gene_name, record.host_organism, record.host_taxonomy_id, record.keywords, record.molecule_type, record.organelle, record.organism, record.organism_classification, record.protein_existence, record.references, record.seqinfo, record.sequence, record.sequence_length, record.sequence_update, record.taxonomy_id]
for i in range(len(info)):
    print(dir(record)[i+26],":" ,end=" ")
    if info[i] != record.references: 
        print(info[i])
        print()
    else:
        for k in range(len(info[i])):
            print()
            print("Authors:" ,end=" ")
            print(info[i][k].authors)
            print("comments:" ,end=" ")
            print(info[i][k].comments)
            print("location:" ,end=" ")
            print(info[i][k].location)
            print("number: ", end=" ")
            print(info[i][k].number)
            print("positions:" ,end=" ")
            print(info[i][k].positions)
            print("references:" ,end=" ")
            print(info[i][k].references)
            print("title:" ,end=" ")
            print(info[i][k].title)
            print()
            


#dir(record)
#%%
#retira os dados do NCBI CDD para essa proteína 
handle=Entrez.efetch(db='protein', rettype='gb', retmode='text', id='P63638')
Glu_racemase=SeqIO.read(handle, 'genbank')
SeqIO.write(Glu_racemase,"Glutamate_racemase.gbk","genbank")
handle.close()
for i in Glu_racemase.features:
    print(i)
    print(i.qualifiers)
    

#%%
#retira os dados das proteínas do Uniprot
#'SA0457', UDP-N-acetylglucosamine 
#https://www.uniprot.org/uniprot/Q7A7B4
# annotation score: 5/5
 
handle = ExPASy.get_sprot_raw('Q7A7B4')
record = SwissProt.read(handle)
#print(record.description)
#for ref in record.references:
#    print("authors:", ref.authors)
#    print("title:", ref.title)
#print(record.organism_classification)

#dir(record)

#%%
info=[record.accessions, record.annotation_update, record.comments, record.created, record.cross_references, record.data_class, record.description, record.entry_name, record.features, record.gene_name, record.host_organism, record.host_taxonomy_id, record.keywords, record.molecule_type, record.organelle, record.organism, record.organism_classification, record.protein_existence, record.references, record.seqinfo, record.sequence, record.sequence_length, record.sequence_update, record.taxonomy_id]
for i in range(len(info)):
    print(dir(record)[i+26],":" ,end=" ")
    if info[i] != record.references: 
        print(info[i])
        print()
    else:
        for k in range(len(info[i])):
            print()
            print("Authors:" ,end=" ")
            print(info[i][k].authors)
            print("comments:" ,end=" ")
            print(info[i][k].comments)
            print("location:" ,end=" ")
            print(info[i][k].location)
            print("number: ", end=" ")
            print(info[i][k].number)
            print("positions:" ,end=" ")
            print(info[i][k].positions)
            print("references:" ,end=" ")
            print(info[i][k].references)
            print("title:" ,end=" ")
            print(info[i][k].title)
            print()
            

#%%
#retira os dados do NCBI CDD para essa proteína 
handle=Entrez.efetch(db='protein', rettype='gb', retmode='text', id='Q7A7B4')
glmU=SeqIO.read(handle, 'genbank')
SeqIO.write(glmU,"BifunctionalProteinGlmU.gbk","genbank")
handle.close()
for i in glmU.features:
    print(i)
    print(i.qualifiers)
