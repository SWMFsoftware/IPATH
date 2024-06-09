       
       integer       ranks,        line,        i,            j
       parameter     ( ranks =4, line =40*50)
       real*8        fp_data (4, line, ranks),
     1               fp_total (4, line)
       
       character(len=1024) :: filename
      
       do i =1, ranks
         j=i-1+10
         write(filename,"(A3,I2.2)") "fp_", i-1
         open(j, file=filename)
       enddo

       do i=1, ranks
         j=i-1+10
         read(j,*) fp_data(1:4, 1:line, i)
         close(j)  
       enddo

       open(21, file="fp_total", form="formatted")
          
       do j = 1, line        
          do i = 1,4
              fp_total(i,j) = sum(fp_data(i, j, 1:ranks))/ranks

          enddo
c          write(21, "(I2, 3ES14.7)") int(fe_total(1,j)), 
          write(21, 1000) int(fp_total(1,j)), 
     1                      fp_total(2:4,j)
       enddo
       
       close(21)
1000   format(I4, ES25.14E3, ES25.14E3, ES25.14E3)
      end
