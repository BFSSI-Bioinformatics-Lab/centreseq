Examples
~~~~~~~~


Basic
^^^^^

Given a directory with the following 10 complete genome assemblies:

::

   /home/user/assemblies/
    ├── Amager.1684.BMH.fasta
    ├── Braenderup.276.BMH.fasta
    ├── Citro.freundii.ID133751.fasta
    ├── Citro.sp.ID135111.fasta
    ├── Citro.werkmanii.ID060152.fasta
    ├── Daytona.100.BMH.fasta
    ├── Duesseldorf.2372.BMH.fasta
    ├── Enteritidis.GC138550.fasta
    ├── Enteritidis.GC148474.fasta
    └── Fresno.2879.BMH.fasta

We can run the following command to run centreseq against our 10 input genomes with 16 cores with all output being
routed to `/home/user/centreseq_output`:

::

    centreseq core -f /home/user/assemblies/ -o /home/user/centreseq_output -n 16

