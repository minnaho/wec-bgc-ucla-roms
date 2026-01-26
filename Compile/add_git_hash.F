      subroutine add_git_hash(ncid) 

      use roms_read_write
      use nc_read_write
      use netcdf

      implicit none
      integer :: ncid,ierr

      git_hash="2a4e9107ffc771d2833fbea3e08f1e169759b5ff"

      ierr=nf90_put_att(ncid, nf90_global, 'git_version',git_hash) 
      if (ierr/=nf90_noerr)
     &    call handle_ierr(ierr,'add_git hash: ')


      end subroutine add_git_hash

