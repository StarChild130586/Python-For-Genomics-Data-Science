#!/usr/bin/env python
# coding: utf-8

# In[1]:


try:
    f=open('sequence.fasta')
except IOError:
    print("the fasta file does not exists")


# In[2]:


pwd


# In[3]:


f


# In[4]:


seqs={}


# In[5]:


for line in f:
 
    line=line.rstrip()
    print(line)
    if line[0]=='>' or  line.startswith('>'):
        words=line.split()
        name=words[0][1:]
        seqs[name]=""
    else:
        seqs[name]=seqs[name]+line
close(f)


# In[6]:


name


# In[7]:


seqs[name]


# In[8]:


for name,seq in seqs.items():
   print(name,seq)


# In[9]:


import Bio


# In[10]:


print(Bio.__version__)


# In[11]:


from  Bio.Blast import  NCBIWWW


# In[12]:


fasta_string='CATGCTACGGTGCTAAAAGCATTACGCCCTATAGTGATTTTCGAGACATACTGTGTTTTTAAATATAGTATTGCC'


# In[13]:


fasta_string


# In[ ]:


len(fasta_string)


# In[ ]:


result_handle=NCBIWWW.qblast("blastn","nt",fasta_string)

				Do	a	BLAST	search	using	the	QBLAST	server	at	NCBI.	
				Supports	all	parameters	of	the	qblast	API	for	Put	and	Get.	
				Some	useful	parameters:	
					-	program								blastn,	blastp,	blastx,	tblastn,	or	tblastx	(lower	case)	
					-	database							Which	database	to	search	against	(e.g.	"nr").	
					-	sequence							The	sequence	to	search.	
					-	ncbi_gi								TRUE/FALSE	whether	to	give	'gi'	identiOier.	
					-	descriptions			Number	of	descriptions	to	show.		Def	500.	
					-	alignments					Number	of	alignments	to	show.		Def	500.	
					-	expect									An	expect	value	cutoff.		Def	10.0.	
					-	matrix_name				Specify	an	alt.	matrix	(PAM30,	PAM70,	BLOSUM80,	BLOSUM45).	
					-	Oilter									"none"	turns	off	Oiltering.		Default	no	Oiltering	
					-	format_type				"HTML",	"Text",	"ASN.1",	or	"XML".		Def.	"XML".	
					-	entrez_query			Entrez	query	to	limit	Blast	search	
					-	hitlist_size			Number	of	hits	to	return.	Default	50	
					-	megablast						TRUE/FALSE	whether	to	use	MEga	BLAST	algorithm	(blastn	only)	
					-	service								plain,	psi,	phi,	rpsblast,	megablast	(lower	case 
# In[ ]:


from Bio.Blast import NCBIXML


# In[ ]:


blast_record = NCBIXML.read(result_handle)


# In[ ]:


len(blast_record.alignments)


# In[ ]:


E_VALUE_THRESH = 0.01


# In[ ]:


for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print('****Alignment****')
            print('sequence:', alignment.title)
            print('length:', alignment.length)
            print('e value:', hsp.expect)
            print(hsp.query)
            print(hsp.match)
            print(hsp.sbjct)

Nucleotide sequences and (reverse) complements
# In[14]:


from Bio.Seq import Seq


# In[15]:


my_seq=Seq(fasta_string)


# In[16]:


my_seq


# In[17]:


my_seq.reverse_complement()


# In[18]:


my_seq.complement()


# In[33]:


my_seq='TGGGCCTCATATTTATCCTATATACCATGTTCGTATGGTGGCGCGATGTTCTACGTGAATCCACGTTCGAAGGACATCATACCAAAGTCGTACAATTAGGACCTCGATATGGTTTTATTCTGTTTATCGTATCGGAGGTTATGTTCTTTTTTGCTCTTTTTCGGGCTTCTTCTCATTCTTCTTTGGCACCTACGGTAGAG'


# In[34]:


my_seq=Seq(my_seq)


# In[41]:


messenger_rna=my_seq.reverse_complement().transcribe()

Sticking with the same example discussed in the transcription section above, now let’s translate this mRNA into the corresponding protein sequence - again taking advantage of one of the Seq object’s biological methods:
# In[42]:


messenger_rna.translate()


# In[43]:


my_seq.translate()


# In[45]:


from Bio.Seq import Seq
my_seq = Seq("AGTACACTGGT")
print(my_seq)
 


# In[46]:


print(my_seq + " - Sequence")
print(my_seq.complement() + " - Complement")
print(my_seq.reverse_complement() + " - Reverse Complement")


# In[47]:


from Bio import SeqIO
count = 0
sequences = [] # Here we are setting up an array to save our sequences for the next step

for seq_record in SeqIO.parse("sequence.fasta", "fasta"):
    if (count < 6):
        sequences.append(seq_record)
        print("Id: " + seq_record.id + " \t " + "Length: " + str("{:,d}".format(len(seq_record))) )
        print(repr(seq_record.seq) + "\n")
        count = count + 1


# In[55]:


chr2L=sequences[0].seq

The Seq object has a .count() method, just like a string. Note that this means that like a Python string, this gives a non-overlapping count
# In[56]:


print("AAAAAA".count("AA"))
print(Seq("AAAA").count("AA"))


# In[58]:


print("Length:\t" + str(len(chr2L)))
print("G Count:\t" + str(chr2L.count("G")))


# In[59]:


"""The GC Content of a DNA sequence is important and relates to how stable the molecule will be. We can calculate it manually like this:"""


# In[60]:


print("GC%:\t\t" + str(100 * float((chr2L.count("G") + chr2L.count("C")) / len(chr2L) ) ))


# In[61]:


from Bio.SeqUtils import GC
print("GC% Package:\t" + str(GC(chr2L)))


# In[62]:


print("GgCcSs%:\t" + str(100 * float((chr2L.count("G") + chr2L.count("g") + chr2L.count("C") + chr2L.count("c") + chr2L.count("S") + chr2L.count("s") ) / len(chr2L) ) ))
print("GC% Package:\t" + str(GC(chr2L)))


# In[65]:


chr2LSHORT = chr2L[0:20]
print("Coding DNA: " + chr2LSHORT)
template_dna = chr2LSHORT.reverse_complement()
print("Template DNA: " + template_dna)


# In[66]:


"""Biology Note: (remember by convention nucleotide sequences are normally read from the 5’ to 3’ direction)

Now let’s transcribe the coding strand into the corresponding mRNA, using the Seq object’s built in transcribe method:"""


# In[67]:


messenger_rna = chr2LSHORT.transcribe()
print("Messenger RNA: " + messenger_rna)


# In[68]:


"""Translation¶
Using a new example, let’s translate this mRNA into the corresponding protein sequence - again taking advantage of one of the Seq object’s biological method"""


# In[70]:


print("Protein Sequence: " + messenger_rna.translate())


# In[ ]:




