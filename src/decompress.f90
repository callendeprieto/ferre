subroutine decompress(influx,outflux)

use share, only: dp,scalef,scaled,npix,                                      &
   npca,nelnpca,nvar,totalnpca,f_access,constant,pcaproject,wpca,meanspca

implicit none

!in/out
real(dp),intent(in)	   :: influx(npix)		!synthetic fluxes
real(dp),intent(out)   :: outflux(totalnpca)		!npca expanded flux vector

!locals
integer				   :: i,ii,jj,kk,offset
real(dp)               :: flux(npix)			!tmp data holder 


flux=influx

if (scaled == 1) flux=flux*scalef
flux = flux - constant

do ii=1,nelnpca
  	offset=0
	if (ii > 1) offset=sum(npca(1:ii-1))
	do i=1,npca(ii)
		kk=i+offset
		outflux(kk)=0.0_dp
		do jj=1,nvar 
			outflux(kk) = outflux(kk) + flux((ii-1)*nvar+jj)*wpca(kk,jj)
		enddo
	enddo
enddo		  	
outflux = outflux + meanspca
	

end subroutine decompress
