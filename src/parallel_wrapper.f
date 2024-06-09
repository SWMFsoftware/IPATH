      subroutine initnodes(myrank,nodes)
      include 'mpif.h'
      
      call MPI_Init(ierr)
      call MPI_COMM_rank(MPI_COMM_WORLD,myrank,ierr)
      call MPI_Comm_size(MPI_COMM_WORLD,nodes,ierr)

      end

      subroutine finishnodes
      include 'mpif.h'
      call MPI_FINALIZE(ierr)
      end

      subroutine getrankstr(rankstr,rank)
      integer rank, node
      character rankstr(*)
      inta = ichar('0')
c      trank = rank + node*8 
      rankstr(1)=char(inta+rank/10)
      rankstr(2)=char(inta+mod(rank,10))
      print*, " === rankstr(1) and rank === ", rankstr(1), rank
      print*, " === rankstr(2) is rank === ", rankstr(2), rank
      end

