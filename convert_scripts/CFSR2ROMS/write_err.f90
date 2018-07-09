subroutine write_err(status)
  integer status
  write(*,*) 'Error, message = ',status
  stop "stopped"
end subroutine write_err
