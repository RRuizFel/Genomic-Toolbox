import collections

class DNATools:
    def __init__(self):
        print('DNA Tools ready for use')
    #Nucleotide Bases
    Nucleotides = ["A", "T", "G", "C",]
    #RNA Codons 
    RNA_codons = {
        'UUU':'F',      'CUU':'L',   'AUU':'I',   'GUU':'V', 'UUC':'F',   'CUC':'L',  'AUC':'I',    'GUC':'V',
    'UUA':'L', 'CUA':'L',  'AUA':'I',      'GUA':'V','UUG':'L',    'CUG':'L',      'AUG':'M',      'GUG':'V',
    'UCU':'S', 'CCU':'P',   'ACU':'T',      'GCU':'A','UCC':'S',      'CCC':'P',      'ACC':'T',      'GCC':'A',
    'UCA':'S',  'CCA':'P',   'ACA':'T',      'GCA':'A','UCG':'S',   'CCG':'P',      'ACG':'T',      'GCG':'A',
    'UAU':'Y',    'CAU':'H',   'AAU':'N',      'GAU':'D','UAC':'Y',      'CAC':'H',      'AAC':'N',      'GAC':'D',
    'UAA':'*', 'CAA':'Q',   'AAA':'K',      'GAA':'E','UAG':'*',   'CAG':'Q',      'AAG':'K',      'GAG':'E',
    'UGU':'C',   'CGU':'R',  'AGU':'S',      'GGU':'G','UGC':'C',      'CGC':'R',      'AGC':'S',      'GGC':'G',
    'UGA':'*',  'CGA':'R',   'AGA':'R',      'GGA':'G','UGG':'W',     'CGG':'R',      'AGG':'R',      'GGG':'G'
    }
    #Monoistopic mass table for amino acids
    Monoistopic_mass_table = {'A':   71.03711, 'C':103.00919, 'D':115.02694, 'E':129.04259,'F':147.06841, 
    'G':57.02146, 'H':137.05891, 'I':113.08406, 'K':128.09496,'L':113.08406, 'M':131.04049, 'N':114.04293, 
    'P':97.05276, 'Q':128.05858,'R':156.10111, 'S':87.03203, 'T':101.04768, 'V':99.06841, 'W':186.07931, 'Y':163.06333 ,
    }

    #Break up content of fafsta files into a dict
    def readFafstaFiles(file):
        with open(file, 'r') as f:
            lines = f.readlines()
        # fafstaDict = {}
        # fafstaLabel = ''
        # for line in lines:
        #     if '>' in line:
        #         fafstaLabel = line
        #         fafstaDict[fafstaLabel] = ''
        #     else:
        #         fafstaDict[fafstaLabel] += line
        # return fafstaDict
                
        
    # ENSURE DNA IS VALID
    def Validate_Seq(dna_Seq):
        tmpSeq = dna_Seq.upper()
        for nuc in tmpSeq:
            if nuc not in DNATools.Nucleotides:
                return False
        return tmpSeq

    # COUNT NUCLEOTIDE FREQUENCY
    def Nucleotide_Freq(dna_seq):
        tmpNucDict = {"A":0, "C":0, "G":0, "T":0}
        for nuc in dna_seq:
            if nuc in tmpNucDict:
                tmpNucDict[nuc] +=1
        return tmpNucDict
        # return dict(collections.Counter(dna_seq))
            #Returns an unfixed order

    # TRANSCRIBE DNA INTO RNA
    # RETURN BASE COUNT
    def DNAtoRNA(dna_seq):
        rnaString = dna_seq.replace('T', 'U')
        print(rnaString)
        tmpNucDict = {"A":0, "C":0, "G":0, "U":0}
        for nuc in rnaString:
            if nuc in tmpNucDict:
                tmpNucDict[nuc] +=1
        return tmpNucDict

    #REVERSE COMPLEMENT OF A DNA STRING
    def reverse_Complement(dna_seq):
        complementDict = {'A':'T','T':'A','G':'C','C':'G'}
        revSeq = dna_seq[::-1]
        revCompStr = ''.join(complementDict.get(nuc,nuc) for nuc in revSeq)
        #revCompStr = ''.join(complementDict.get(nuc,nuc) for nuc in reverse(dna_seq))
        return revCompStr

    #GC CONTENT NON-PERCENT CALCULATOR
    def GC_content(seq):
        G_content = seq.count('G')
        C_content = seq.count('C')
        results = (G_content + C_content)/len(seq)*100
        return round(results, 4)

    #GC CONTENT PERCENT CALCULATOR
    def GC_content_percent(seq):
        G_content = seq.count('G')
        C_content = seq.count('C')
        results = ((G_content + C_content)/len(seq))*100
        GC_CalcPerc = '{:.2f}%'.format(results)
        return GC_CalcPerc

    # HIGHEST GC CONTENT OF VARIOUS FAFSTA DNA SEQUENCES
    # GC CONTENT = REVERSE COMPLEMENT GC CONTENT
    def Max_GC_Content(fafstaFile):
        with open(fafstaFile, 'r') as f:
            lines = f.readlines()
        fafstaDict = {}
        fafstaLabel = ''
        for line in lines:
            if '>' in line:
                fafstaLabel = line
                fafstaDict[fafstaLabel] = ''
            else:
                fafstaDict[fafstaLabel] += line
        resultDict = {key: DNATools.GC_content_percent(value) for (key,value) in fafstaDict.items()}
        Max_GC_Key = max(resultDict, key = resultDict.get)
        #From the resultDict, iterate through every key and return the key with max value
        #Key to call value from resultDict later on
        print(f'{Max_GC_Key}{resultDict[Max_GC_Key]}')


    # TRANSLATE RNA INTO PROTEIN
    def RNAtoProtein(rnaSeq):
        protein = ''
        for i in range(0, len(rnaSeq), 3):
            if DNATools.RNA_codons[rnaSeq[i:i+3]] != '*': 
                protein += DNATools.RNA_codons[rnaSeq[i:i+3]]
            else:
                return protein
        #SOLUTION 2                 Goes through entire sequence regardless of stop codon
        # for i in range(0, len(rnaSeq), 3):
        #     codon = rnaSeq[i:i+3]
        #     if codon in RNA_codons:
        #         protein += RNA_codons[codon]
        # return protein
            

    #FIND MOTIFS IN DNA
    def DNAmotifs(DNAseq):
        motif = input('What motif are you looking for in this DNA sequence?')
        positions = []
        for i in range(len(DNAseq)):
            if DNAseq[i] == motif[0]:
                if DNAseq[i : i + len(motif)] == motif:
                    positions.append(i+1)
        return positions 
        # print(*positions)
        
    #PROTEIN MASS OF A PROTEIN
    def ProteinMass(ProteinSeq):
        mass = 0.0
        for base in ProteinSeq:
            mass += DNATools.Monoistopic_mass_table.get(base, 0.0)
        massRounded = round(mass, 3)
        return massRounded
        
    # POINT MUTATION COUNTER, DETECTS DIFFERENT BASES IN 2 SEQUENCES
    def HammingDistance():
        counter = 0
        seq1 = input('Input the first sequence: ')
        seq2 = input('Input the second sequence: ')
        if len(seq1) != len(seq2):
            UserResponse = input('Sequences are of different lengths\n To proceed type "Yes", otherwise "No": ').lower()
            if UserResponse == 'yes':
                DNATools.HammingDistance_DifferentSizes(seq1, seq2)
            elif UserResponse == 'no':
                return None
            else:
                UserResponse = input('Invalid input\n To proceed type "Yes", otherwise "No": ').lower()
        for char1, char2 in zip(seq1, seq2):
            if char1 != char2:
                counter += 1
        return counter         
    def HammingDistance_DifferentSizes(seq1, seq2):
        counter = 0
        maxSeq = max(len(seq1), len(seq2))
        for i in range(maxSeq):
            if i >= len(seq1) or i >= len(seq2) or seq1[i] != seq2[i]:
                counter += 1
        return counter

    #LOCATION OF DIFFERENT BASES
    def mutationLocation():
        positions = []
        seq1 = input('Input the first sequence: ')
        seq2 = input('Input the second sequence: ')
        if len(seq1) != len(seq2):
            UserResponse = input('Sequences are of different lengths\n To proceed type "Yes", otherwise "No": ').lower()
            if UserResponse == 'yes':
                DNATools.muationLocation_DifferentSizes(seq1, seq2)
            elif UserResponse == 'no':
                return None
            else:
                UserResponse = input('Invalid input\n To proceed type "Yes", otherwise "No": ').lower()
        for i in range(min(len(seq1), len(seq2))):
            if seq1[i] != seq2[i]:
                positions.append(i+1)
        return positions
    def muationLocation_DifferentSizes(seq1, seq2):
        positions = []
        for i, (base1, base2) in enumerate(zip(seq1, seq2)):
            if base1 != base2:
                positions.append(i + 1)
        max_length = max(len(seq1), len(seq2))
        if len(seq1) < max_length:
            for j in range(len(seq1), max_length):
                positions.append(j)
        elif len(seq2) < max_length:
            for j in range(len(seq2), max_length):
                positions.append(j)
        return positions
        
        
    # #FIBONACCI PAIR POPULATION
    def Fibonacci_pair_population(num_months, num_Offspring):
        Months_Offspring = [0, 1]
        if num_months == 1:
            return 1
        elif num_months == 2:
            return num_Offspring
        for i in range(2, num_months + 1):
            Months_Offspring.append(Months_Offspring[i-1] + ((Months_Offspring[i - 2]) * num_Offspring))
        # BELOW returns entire Fib loop
        #return Months_Offspring
        return Months_Offspring[-1]

    # # MENDELIAN INHERITANCE
    # def mendels_first_law(k, m, n):
    #     total_pop = k + m + n

    #Returns the two FAFSTA files that have overlapping bases 
    # def OverlappGraphs(FafstaFile):
    #     f