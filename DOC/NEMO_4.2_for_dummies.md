# NEMO4.2 for dummies
## Introduction:
NEMO4.2 is a major release of NEMO. Compared to NEMO 4.0.x, there are a lot of new features that have a deep impact on the model architecture 
and on the way to use it for a simulation production.  This document was written when discovering the new features and tries to explain 
the changes from the user point of view in a pragmatic way.  Each specific point is treated in a chapter.  Chapters are written when
digging into a particular topics. The order of the chapters is not necessarily important.

## Use of macro for DO loop:
### The main ideas:
This is probably the most visible novelty when you open a 4.2 fortran module, compared to 4.0 or older release. At compile time, all the CPP
macros are replace by the due imbedded fortran loops, using definitions given in `do_loop_substitute.h90`, which is reproduced *in extenso* below:

```fortran
#if defined show_comments
! These comments are not intended to be retained during preprocessing; i.e. do not define "show_comments"
!!----------------------------------------------------------------------
!! NEMO/OCE 4.x , NEMO Consortium (2020)
!! Software governed by the CeCILL license (see ./LICENSE)
!!----------------------------------------------------------------------
! This header file contains preprocessor definitions and macros used in the do-loop substitutions introduced
! between version 4.0 and 4.2. The primary aim of these macros is to assist in future applications of tiling
! to improve performance. This is expected to be achieved by alternative versions of these macros in selected
! locations. The initial introduction of these macros simply replaced all identifiable nested 2D- and 3D-loops
! with single line statements (and adjusts indenting accordingly). Do loops were identifiable if they comformed
! to either:
!                                       DO jk = ....
!   DO jj = ....                           DO jj = ...
!      DO ji = ....                           DO ji = ...
!         .                   OR                 .
!         .                                      .
!     END DO                                  END DO
!   END DO                                 END DO
!                                       END DO
! and white-space variants thereof.
!
! Additionally, only loops with recognised jj and ji loops limits were treated; these were:
! Lower limits of 1, 2 or fs_2
! Upper limits of jpi, jpim1 or fs_jpim1 (for ji) or jpj, jpjm1 or fs_jpjm1 (for jj)
!
! The macro naming convention takes the form: DO_2D( L, R, B, T) where:
!   L is the Left   offset from the PE's inner domain;
!   R is the Right  offset from the PE's inner domain
!   B is the Bottom offset from the PE's inner domain;
!   T is the Top    offset from the PE's inner domain;
!
! So, given an inner domain of 2,jpim1 and 2,jpjm1, a typical example would replace:
!
!   DO jj = 2, jpj
!      DO ji = 1, jpim1
!         .
!         .
!      END DO
!   END DO
!
! with:
!
!   DO_2D( 1, 0, 0, 1 )
!      .
!      .
!   END_2D
!
! similar conventions apply to the 3D loops macros. jk loop limits are retained through macro arguments
! and are not restricted. This includes the possibility of strides for which an extra set of DO_3DS
! macros are defined.
!
! In the following definitions the inner PE domain is defined by start indices of (Nis0, Njs0) and end
! indices of (Nie0, Nje0) where:
!
! Nis0 =   1 + nn_hls     Njs0 =   1 + nn_hls
! Nie0 = jpi - nn_hls     Nje0 = jpj - nn_hls
!
#endif

#define DO_2D(L, R, B, T) DO jj = ntsj-(B), ntej+(T) ; DO ji = ntsi-(L), ntei+(R)
#define DO_2D_OVR(L, R, B, T) DO_2D(L-(L+R)*nthl, R-(R+L)*nthr, B-(B+T)*nthb, T-(T+B)*ntht)
#define A1Di(H) ntsi-(H):ntei+(H)
#define A1Dj(H) ntsj-(H):ntej+(H)
#define A2D(H) A1Di(H),A1Dj(H)
#define A1Di_T(T) (ntsi-nn_hls-1)*T+1:
#define A1Dj_T(T) (ntsj-nn_hls-1)*T+1:
#define A2D_T(T) A1Di_T(T),A1Dj_T(T)
#define JPK  :
#define JPTS  :
#define KJPT  :

#define DO_3D(L, R, B, T, ks, ke) DO jk = ks, ke ; DO_2D(L, R, B, T)
#define DO_3D_OVR(L, R, B, T, ks, ke) DO jk = ks, ke ; DO_2D_OVR(L, R, B, T)

#define DO_3DS(L, R, B, T, ks, ke, ki) DO jk = ks, ke, ki ; DO_2D(L, R, B, T)
#define DO_3DS_OVR(L, R, B, T, ks, ke, ki) DO jk = ks, ke, ki ; DO_2D_OVR(L, R, B, T)

#define END_2D   END DO   ;   END DO
#define END_3D   END DO   ;   END DO   ;   END DO
```

  * The key points to retain are the following:

    ```fortran
! The macro naming convention takes the form: DO_2D( L, R, B, T) where:
!   L is the Left   offset from the PE's inner domain;
!   R is the Right  offset from the PE's inner domain
!   B is the Bottom offset from the PE's inner domain;
!   T is the Top    offset from the PE's inner domain;
    ```

    and a fairly easy to understand example when remembering previous NEMO versions:

    ```
! So, given an inner domain of 2,jpim1 and 2,jpjm1, a typical example would replace:
!
!   DO jj = 2, jpj
!      DO ji = 1, jpim1
!         .
!         .
!      END DO
!   END DO
!
! with:
!
!   DO_2D( 1, 0, 0, 1 )
!      .
!      .
!   END_2D

    ```
#### Definition of the inner domain:
  This is the central part of the DO loop macro definition: 

## CPP keys : actual list and purposes:

## namelists_ref

## XIOS fields

## Tiling : ? 
