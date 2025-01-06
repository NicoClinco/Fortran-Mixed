subroutine template()
  use OAD_tape
  use OAD_rev

  !The template routine used for transforming
  !"getLayerIndexes" subroutine

  !$TEMPLATE_PRAGMA_DECLARATIONS

  integer iaddr
  external iaddr

  if (our_rev_mode%plain) then
     gindxs(1) = (ilayer-1)*nin*nout + (i-1)*nin+j
     gindxs(2) = (nlayers)*nin*nout+(ilayer-1)*nout + i
  end if
  if (our_rev_mode%tape) then
     gindxs(1) = (ilayer-1)*nin*nout + (i-1)*nin+j
     gindxs(2) = (nlayers)*nin*nout+(ilayer-1)*nout + i
  end if
  if (our_rev_mode%adjoint) then
     gindxs(1) = (ilayer-1)*nin*nout + (i-1)*nin+j
     gindxs(2) = (nlayers)*nin*nout+(ilayer-1)*nout + i
  end if
end subroutine template
