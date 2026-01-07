from libc.stdio cimport FILE

cdef extern from "aragorn.h" nogil:

    const size_t STRLEN
    const size_t STRLENM1
    const size_t SHORTSTRLEN
    const size_t SHORTSTRLENM1

    const size_t NS

    const size_t MAXGCMOD
    const size_t MAMMAL_MT
    const size_t NGENECODE

    const size_t MININTRONLEN
    const size_t MAXINTRONLEN
    const size_t MINCTRNALEN
    const size_t MAXCTRNALEN
    const size_t MINTRNALEN
    const size_t MAXTRNALEN
    const size_t MAXETRNALEN

    const size_t MINTAGDIST
    const size_t MAXTAGDIST

    const size_t MAXTMRNALEN
    const size_t TSWEEP
    const size_t WRAP
    const size_t NPTAG
    const size_t MAXAGENELEN

    const size_t NA
    const size_t ND
    const size_t NT
    const size_t NH
    const size_t NTH
    const size_t NC
    const size_t NGFT
    const size_t NTAG
    const size_t NTAGMAX
    const size_t LSEQ
    const size_t ATBOND

    cdef bint space(char c)
    cdef long int sq(long int pos)

    cdef enum base:
        INSERT
        TERM
        Adenine
        Cytosine
        Guanine
        Thymine
        AMBIG
        NOBASE

    cdef enum data_type:
        FASTA
        GENBANK

    cdef enum gencode:
        METAZOAN_MT
        STANDARD
        VERTEBRATE_MT

    cdef enum gene_type:
        noGENE
        tRNA
        tmRNA
        srpRNA
        rRNA
        CDS

    ctypedef struct annotated_gene:
        long int start
        long int stop
        int comp
        long int antistart
        long int antistop
        gene_type genetype
        bint pseudogene
        bint permuted
        bint detected
        char species[SHORTSTRLEN]

    ctypedef struct data_set:
        char filename[80]
        FILE *f
        char seqname[STRLEN]
        bint bugmode
        data_type datatype
        double gc
        long int filepointer
        long int ps
        long int psmax
        long int seqstart
        long int seqstartoff
        long int nextseq
        long int nextseqoff
        int ns
        int nf
        long int aseqlen
        int nagene[NS]
        annotated_gene gene[NGFT]

    ctypedef struct gene:
        char name[100]
        int seq[MAXTRNALEN + 1]
        int eseq[MAXETRNALEN + 1]
        int *ps
        int nbase
        int comp
        long int start
        long int stop
        int astem1
        int astem2
        int aatail
        int spacer1
        int spacer2
        int dstem
        int dloop
        int cstem
        int cloop
        int intron
        int nintron
        int anticodon
        int var
        int varbp
        int tstem
        int tloop
        int genetype
        double energy
        int asst
        int tps
        int tpe
        int annotation
        int annosc

    ctypedef struct csw:
        char genetypename[NS][10]
        FILE *f
        int batch
        int batchfullspecies
        int repeatsn
        int trna
        int tmrna
        int srprna
        int cds
        int mtrna
        int tvloop
        int cloop7
        int peptide
        int geneticcode
        int ngcmod
        int gcmod[MAXGCMOD]
        int gcfix
        int discrim
        int extastem
        int tarm
        int tagthresh
        int tarmlength
        int showconfig
        int libflag
        bint verbose
        int linear
        int both
        int reportpseudogenes
        int energydisp
        int secstructdisp
        int seqdisp
        int aataildisp
        int aataildiv
        int sp1max
        int sp2min
        int sp2max
        int mtxdetect
        int mtcdsscan
        int mtcompov
        int matchacceptor
        int maxintronlen
        int minintronlen
        int minintronlenreport
        int ioverlay
        int ifixedpos
        int ireportminintronlen
        int tmstrict
        int iamismatch
        int loffset
        int roffset
        long int start
        int comp
        gene* genes
        int genespace
        int srpspace
        int ngene[NS]
        int nps
        int annotated
        int dispmatch
        int updatetmrnatags
        int tagend
        int trnalenmisthresh
        int tmrnalenmisthresh
        int nagene[NS]
        int nafn[NS]
        int nafp[NS]
        int natfpd
        int natfptv
        int lacds
        int ldcds
        long int nabase
        double reportpsthresh
        double threshlevel
        double trnathresh
        double ttscanthresh
        double ttarmthresh
        double tdarmthresh
        double tastemthresh
        double tascanthresh
        double mttthresh
        double mtdthresh
        double mtdtthresh
        double mttarmthresh
        double mtdarmthresh
        double tmrnathresh
        double tmathresh
        double tmcthresh
        double tmcathresh
        double tmrthresh
        double srpthresh
        double cdsthresh
        double eref[NS]
        int tmrna_struct[200]

    int[3][6][6] lbp
    int[6][6] bp
    int[6][6] wbp
    int[6][6] wcbp
    int[6][6] gc
    int[6][6] gt
    int[6][6] at
    int[6][6] tt
    int[6][6] stemterm
    int[6][6] aastemterm
    int[6][6] ggstemterm
    int[6][6] assymst
    int[6][6] assymat
    int[6][6] stackbp
    int[6][6] ggstackbp
    int[6][6] ggbp
    int[6][6] gabp
    int[6][6] assymagbp
    int[6][6] stembp
    int[6][6] ggstembp
    int[6][6] gastembp
    int[6][6] vbp
    int[256] map
    # int[mtNTM][4] tandemid
    # double[mtNTM] tandem_em
    # double[6][6] send_em
    # double[6][6] ssend_em
    # int[6][6] neighbour_map
    # double[2][6][6] neighbour_em
    # unsigned int[6][6] btmap
    # double[6][6] bem
    # int[3][64][6] mt_discrim
    # char[NAMINOACID + 1] aapolarity
    # char[NAMINOACID + 1] aaletter
    # char[NAMINOACID][20] aaname
    # char[4] ambig_aaname
    # int[NGENECODE][64] aamap
    # char[NHELPLINE][81] helpmenu
    # tmrna_tag_entry[NTAGMAX] tagdatabase

    char* aa(int* anticodon, csw* sw)
    void change_thresholds(csw *sw, double psthresh)
    void bopt_fastafile(data_set *d, csw *sw)
    void batch_gene_set(data_set* d, int nt, csw* sw)
    char cbase(int c)
    char cpbase(int c)
    int gene_sort(data_set *d, int nt, int* sort, csw *sw)
    void init_gene(gene* ts, int nstart, int nstop)
    char ltranslate(int *codon, gene *t, csw *sw)
    int move_forward(data_set *d)
    double nenergy(gene *t, csw *sw)
    void remove_overlapping_trna(data_set *d, int nt, csw *sw)
    void report_new_tmrna_tags(csw *sw)
    void sense_switch(int *seq1, int *seq2, int lseq)
    int seq_init(data_set *d, csw *sw)
    int seqlen(gene* t)
    int tmioptimise(data_set *d, int *seq, int lseq, int nts, csw *sw)
    void update_tmrna_tag_database(gene* ts, int nt, csw *sw)