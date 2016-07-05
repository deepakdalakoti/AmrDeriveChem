#include "CONSTANTS.H"
#include "convertDRM19toC3H8_F.H"

      subroutine FORT_CONVERT_DRM19_TO_C3H8 (t, y)
      implicit none
      double precision t, y(*)
      integer N_old,Nvals_new,Nspec_new,i,Nvals_newMAX
      parameter (Nvals_newMAX = 50)
      double precision x, frac, phr_old, phr_new, sum
      double precision x_data_old(251), tmp(Nvals_newMAX)
      double precision t_data_old(251)
      logical debug
      data N_old /251/
      data phr_old /1.011568d0/
      data phr_new /0.637726d0/ ! from .7 Q profile
      !
      ! X data from drm19_080.f
      !
      data x_data_old(  1) /   0.000000000000D+00/
      data x_data_old(  2) /   0.500000000000D+00/
      data x_data_old(  3) /   0.731625000000D+00/
      data x_data_old(  4) /   0.847437500000D+00/
      data x_data_old(  5) /   0.876390625000D+00/
      data x_data_old(  6) /   0.890867187500D+00/
      data x_data_old(  7) /   0.905343750000D+00/
      data x_data_old(  8) /   0.919820312500D+00/
      data x_data_old(  9) /   0.927058593750D+00/
      data x_data_old( 10) /   0.934296875000D+00/
      data x_data_old( 11) /   0.937916015625D+00/
      data x_data_old( 12) /   0.941535156250D+00/
      data x_data_old( 13) /   0.945154296875D+00/
      data x_data_old( 14) /   0.948773437500D+00/
      data x_data_old( 15) /   0.952392578125D+00/
      data x_data_old( 16) /   0.954202148437D+00/
      data x_data_old( 17) /   0.956011718750D+00/
      data x_data_old( 18) /   0.957821289062D+00/
      data x_data_old( 19) /   0.959630859375D+00/
      data x_data_old( 20) /   0.961440429687D+00/
      data x_data_old( 21) /   0.963250000000D+00/
      data x_data_old( 22) /   0.963537109375D+00/
      data x_data_old( 23) /   0.963824218750D+00/
      data x_data_old( 24) /   0.964398437500D+00/
      data x_data_old( 25) /   0.964972656250D+00/
      data x_data_old( 26) /   0.965546875000D+00/
      data x_data_old( 27) /   0.966695312500D+00/
      data x_data_old( 28) /   0.967843750000D+00/
      data x_data_old( 29) /   0.968417968750D+00/
      data x_data_old( 30) /   0.968992187500D+00/
      data x_data_old( 31) /   0.969566406250D+00/
      data x_data_old( 32) /   0.970140625000D+00/
      data x_data_old( 33) /   0.970714843750D+00/
      data x_data_old( 34) /   0.971289062500D+00/
      data x_data_old( 35) /   0.971863281250D+00/
      data x_data_old( 36) /   0.972437500000D+00/
      data x_data_old( 37) /   0.973011718750D+00/
      data x_data_old( 38) /   0.973585937500D+00/
      data x_data_old( 39) /   0.974160156250D+00/
      data x_data_old( 40) /   0.974734375000D+00/
      data x_data_old( 41) /   0.975308593750D+00/
      data x_data_old( 42) /   0.975882812500D+00/
      data x_data_old( 43) /   0.976457031250D+00/
      data x_data_old( 44) /   0.977031250000D+00/
      data x_data_old( 45) /   0.977605468750D+00/
      data x_data_old( 46) /   0.978179687500D+00/
      data x_data_old( 47) /   0.978466796875D+00/
      data x_data_old( 48) /   0.978753906250D+00/
      data x_data_old( 49) /   0.979041015625D+00/
      data x_data_old( 50) /   0.979328125000D+00/
      data x_data_old( 51) /   0.979615234375D+00/
      data x_data_old( 52) /   0.979902343750D+00/
      data x_data_old( 53) /   0.980189453125D+00/
      data x_data_old( 54) /   0.980476562500D+00/
      data x_data_old( 55) /   0.980763671875D+00/
      data x_data_old( 56) /   0.981050781250D+00/
      data x_data_old( 57) /   0.981337890625D+00/
      data x_data_old( 58) /   0.981625000000D+00/
      data x_data_old( 59) /   0.981768554688D+00/
      data x_data_old( 60) /   0.981912109375D+00/
      data x_data_old( 61) /   0.982055664062D+00/
      data x_data_old( 62) /   0.982199218750D+00/
      data x_data_old( 63) /   0.982486328125D+00/
      data x_data_old( 64) /   0.982773437500D+00/
      data x_data_old( 65) /   0.983060546875D+00/
      data x_data_old( 66) /   0.983347656250D+00/
      data x_data_old( 67) /   0.983634765625D+00/
      data x_data_old( 68) /   0.983921875000D+00/
      data x_data_old( 69) /   0.984208984375D+00/
      data x_data_old( 70) /   0.984496093750D+00/
      data x_data_old( 71) /   0.984783203125D+00/
      data x_data_old( 72) /   0.985070312500D+00/
      data x_data_old( 73) /   0.985357421875D+00/
      data x_data_old( 74) /   0.985644531250D+00/
      data x_data_old( 75) /   0.985931640625D+00/
      data x_data_old( 76) /   0.986218750000D+00/
      data x_data_old( 77) /   0.986792968750D+00/
      data x_data_old( 78) /   0.987367187500D+00/
      data x_data_old( 79) /   0.987941406250D+00/
      data x_data_old( 80) /   0.988515625000D+00/
      data x_data_old( 81) /   0.989089843750D+00/
      data x_data_old( 82) /   0.989664062500D+00/
      data x_data_old( 83) /   0.990238281250D+00/
      data x_data_old( 84) /   0.990812500000D+00/
      data x_data_old( 85) /   0.991386718750D+00/
      data x_data_old( 86) /   0.991673828125D+00/
      data x_data_old( 87) /   0.991960937500D+00/
      data x_data_old( 88) /   0.992535156250D+00/
      data x_data_old( 89) /   0.992822265625D+00/
      data x_data_old( 90) /   0.993109375000D+00/
      data x_data_old( 91) /   0.993396484375D+00/
      data x_data_old( 92) /   0.993683593750D+00/
      data x_data_old( 93) /   0.993970703125D+00/
      data x_data_old( 94) /   0.994257812500D+00/
      data x_data_old( 95) /   0.994544921875D+00/
      data x_data_old( 96) /   0.994832031250D+00/
      data x_data_old( 97) /   0.995119140625D+00/
      data x_data_old( 98) /   0.995406250000D+00/
      data x_data_old( 99) /   0.995693359375D+00/
      data x_data_old(100) /   0.995980468750D+00/
      data x_data_old(101) /   0.996267578125D+00/
      data x_data_old(102) /   0.996554687500D+00/
      data x_data_old(103) /   0.996841796875D+00/
      data x_data_old(104) /   0.997128906250D+00/
      data x_data_old(105) /   0.997416015625D+00/
      data x_data_old(106) /   0.997703125000D+00/
      data x_data_old(107) /   0.997990234375D+00/
      data x_data_old(108) /   0.998277343750D+00/
      data x_data_old(109) /   0.998564453125D+00/
      data x_data_old(110) /   0.998851562500D+00/
      data x_data_old(111) /   0.999138671875D+00/
      data x_data_old(112) /   0.999425781250D+00/
      data x_data_old(113) /   0.999712890625D+00/
      data x_data_old(114) /   0.100000000000D+01/
      data x_data_old(115) /   0.100024414062D+01/
      data x_data_old(116) /   0.100048828125D+01/
      data x_data_old(117) /   0.100073242188D+01/
      data x_data_old(118) /   0.100097656250D+01/
      data x_data_old(119) /   0.100122070312D+01/
      data x_data_old(120) /   0.100146484375D+01/
      data x_data_old(121) /   0.100170898438D+01/
      data x_data_old(122) /   0.100195312500D+01/
      data x_data_old(123) /   0.100219726562D+01/
      data x_data_old(124) /   0.100244140625D+01/
      data x_data_old(125) /   0.100268554688D+01/
      data x_data_old(126) /   0.100292968750D+01/
      data x_data_old(127) /   0.100317382812D+01/
      data x_data_old(128) /   0.100341796875D+01/
      data x_data_old(129) /   0.100366210938D+01/
      data x_data_old(130) /   0.100390625000D+01/
      data x_data_old(131) /   0.100415039062D+01/
      data x_data_old(132) /   0.100439453125D+01/
      data x_data_old(133) /   0.100463867188D+01/
      data x_data_old(134) /   0.100488281250D+01/
      data x_data_old(135) /   0.100512695312D+01/
      data x_data_old(136) /   0.100537109375D+01/
      data x_data_old(137) /   0.100561523438D+01/
      data x_data_old(138) /   0.100585937500D+01/
      data x_data_old(139) /   0.100610351562D+01/
      data x_data_old(140) /   0.100634765625D+01/
      data x_data_old(141) /   0.100659179688D+01/
      data x_data_old(142) /   0.100683593750D+01/
      data x_data_old(143) /   0.100708007812D+01/
      data x_data_old(144) /   0.100732421875D+01/
      data x_data_old(145) /   0.100756835938D+01/
      data x_data_old(146) /   0.100781250000D+01/
      data x_data_old(147) /   0.100805664062D+01/
      data x_data_old(148) /   0.100830078125D+01/
      data x_data_old(149) /   0.100854492188D+01/
      data x_data_old(150) /   0.100878906250D+01/
      data x_data_old(151) /   0.100903320312D+01/
      data x_data_old(152) /   0.100927734375D+01/
      data x_data_old(153) /   0.100952148438D+01/
      data x_data_old(154) /   0.100976562500D+01/
      data x_data_old(155) /   0.101000976562D+01/
      data x_data_old(156) /   0.101025390625D+01/
      data x_data_old(157) /   0.101049804688D+01/
      data x_data_old(158) /   0.101074218750D+01/
      data x_data_old(159) /   0.101098632812D+01/
      data x_data_old(160) /   0.101123046875D+01/
      data x_data_old(161) /   0.101147460938D+01/
      data x_data_old(162) /   0.101171875000D+01/
      data x_data_old(163) /   0.101196289062D+01/
      data x_data_old(164) /   0.101220703125D+01/
      data x_data_old(165) /   0.101245117188D+01/
      data x_data_old(166) /   0.101269531250D+01/
      data x_data_old(167) /   0.101293945312D+01/
      data x_data_old(168) /   0.101318359375D+01/
      data x_data_old(169) /   0.101342773438D+01/
      data x_data_old(170) /   0.101367187500D+01/
      data x_data_old(171) /   0.101391601562D+01/
      data x_data_old(172) /   0.101416015625D+01/
      data x_data_old(173) /   0.101464843750D+01/
      data x_data_old(174) /   0.101513671875D+01/
      data x_data_old(175) /   0.101562500000D+01/
      data x_data_old(176) /   0.101611328125D+01/
      data x_data_old(177) /   0.101660156250D+01/
      data x_data_old(178) /   0.101708984375D+01/
      data x_data_old(179) /   0.101757812500D+01/
      data x_data_old(180) /   0.101806640625D+01/
      data x_data_old(181) /   0.101855468750D+01/
      data x_data_old(182) /   0.101904296875D+01/
      data x_data_old(183) /   0.101953125000D+01/
      data x_data_old(184) /   0.102001953125D+01/
      data x_data_old(185) /   0.102050781250D+01/
      data x_data_old(186) /   0.102099609375D+01/
      data x_data_old(187) /   0.102148437500D+01/
      data x_data_old(188) /   0.102197265625D+01/
      data x_data_old(189) /   0.102246093750D+01/
      data x_data_old(190) /   0.102294921875D+01/
      data x_data_old(191) /   0.102343750000D+01/
      data x_data_old(192) /   0.102392578125D+01/
      data x_data_old(193) /   0.102441406250D+01/
      data x_data_old(194) /   0.102490234375D+01/
      data x_data_old(195) /   0.102539062500D+01/
      data x_data_old(196) /   0.102587890625D+01/
      data x_data_old(197) /   0.102636718750D+01/
      data x_data_old(198) /   0.102685546875D+01/
      data x_data_old(199) /   0.102734375000D+01/
      data x_data_old(200) /   0.102783203125D+01/
      data x_data_old(201) /   0.102832031250D+01/
      data x_data_old(202) /   0.102929687500D+01/
      data x_data_old(203) /   0.103027343750D+01/
      data x_data_old(204) /   0.103125000000D+01/
      data x_data_old(205) /   0.103222656250D+01/
      data x_data_old(206) /   0.103320312500D+01/
      data x_data_old(207) /   0.103417968750D+01/
      data x_data_old(208) /   0.103515625000D+01/
      data x_data_old(209) /   0.103613281250D+01/
      data x_data_old(210) /   0.103710937500D+01/
      data x_data_old(211) /   0.103808593750D+01/
      data x_data_old(212) /   0.103906250000D+01/
      data x_data_old(213) /   0.104101562500D+01/
      data x_data_old(214) /   0.104296875000D+01/
      data x_data_old(215) /   0.104492187500D+01/
      data x_data_old(216) /   0.104687500000D+01/
      data x_data_old(217) /   0.105078125000D+01/
      data x_data_old(218) /   0.105468750000D+01/
      data x_data_old(219) /   0.105859375000D+01/
      data x_data_old(220) /   0.106250000000D+01/
      data x_data_old(221) /   0.106640625000D+01/
      data x_data_old(222) /   0.107031250000D+01/
      data x_data_old(223) /   0.107421875000D+01/
      data x_data_old(224) /   0.107812500000D+01/
      data x_data_old(225) /   0.108203125000D+01/
      data x_data_old(226) /   0.108593750000D+01/
      data x_data_old(227) /   0.108984375000D+01/
      data x_data_old(228) /   0.109375000000D+01/
      data x_data_old(229) /   0.110156250000D+01/
      data x_data_old(230) /   0.110937500000D+01/
      data x_data_old(231) /   0.111718750000D+01/
      data x_data_old(232) /   0.112500000000D+01/
      data x_data_old(233) /   0.113281250000D+01/
      data x_data_old(234) /   0.114062500000D+01/
      data x_data_old(235) /   0.115625000000D+01/
      data x_data_old(236) /   0.117187500000D+01/
      data x_data_old(237) /   0.118750000000D+01/
      data x_data_old(238) /   0.120312500000D+01/
      data x_data_old(239) /   0.121875000000D+01/
      data x_data_old(240) /   0.125000000000D+01/
      data x_data_old(241) /   0.128125000000D+01/
      data x_data_old(242) /   0.131250000000D+01/
      data x_data_old(243) /   0.134375000000D+01/
      data x_data_old(244) /   0.137500000000D+01/
      data x_data_old(245) /   0.143750000000D+01/
      data x_data_old(246) /   0.150000000000D+01/
      data x_data_old(247) /   0.156250000000D+01/
      data x_data_old(248) /   0.162500000000D+01/
      data x_data_old(249) /   0.175000000000D+01/
      data x_data_old(250) /   0.187500000000D+01/
      data x_data_old(251) /   0.200000000000D+01/
      !
      ! T data from drm19_080.f
      !
      data t_data_old(  1) /   0.300000000000D+03/
      data t_data_old(  2) /   0.300000100662D+03/
      data t_data_old(  3) /   0.300002413612D+03/
      data t_data_old(  4) /   0.300030244753D+03/
      data t_data_old(  5) /   0.300104006185D+03/
      data t_data_old(  6) /   0.300247047239D+03/
      data t_data_old(  7) /   0.300664120597D+03/
      data t_data_old(  8) /   0.301876970823D+03/
      data t_data_old(  9) /   0.303346056746D+03/
      data t_data_old( 10) /   0.306194643012D+03/
      data t_data_old( 11) /   0.308610905450D+03/
      data t_data_old( 12) /   0.312132663367D+03/
      data t_data_old( 13) /   0.317235493653D+03/
      data t_data_old( 14) /   0.324568965690D+03/
      data t_data_old( 15) /   0.334991642311D+03/
      data t_data_old( 16) /   0.341793274852D+03/
      data t_data_old( 17) /   0.349928137356D+03/
      data t_data_old( 18) /   0.359599962527D+03/
      data t_data_old( 19) /   0.371023241604D+03/
      data t_data_old( 20) /   0.384417445547D+03/
      data t_data_old( 21) /   0.400000000000D+03/
      data t_data_old( 22) /   0.402699121576D+03/
      data t_data_old( 23) /   0.405464467029D+03/
      data t_data_old( 24) /   0.411195396926D+03/
      data t_data_old( 25) /   0.417197809167D+03/
      data t_data_old( 26) /   0.423478319479D+03/
      data t_data_old( 27) /   0.436881704727D+03/
      data t_data_old( 28) /   0.451434957717D+03/
      data t_data_old( 29) /   0.459163550338D+03/
      data t_data_old( 30) /   0.467202056445D+03/
      data t_data_old( 31) /   0.475554574335D+03/
      data t_data_old( 32) /   0.484224784181D+03/
      data t_data_old( 33) /   0.493215942159D+03/
      data t_data_old( 34) /   0.502530877125D+03/
      data t_data_old( 35) /   0.512171989659D+03/
      data t_data_old( 36) /   0.522141253247D+03/
      data t_data_old( 37) /   0.532440217353D+03/
      data t_data_old( 38) /   0.543070012126D+03/
      data t_data_old( 39) /   0.554031354471D+03/
      data t_data_old( 40) /   0.565324555234D+03/
      data t_data_old( 41) /   0.576949527219D+03/
      data t_data_old( 42) /   0.588905793828D+03/
      data t_data_old( 43) /   0.601192498064D+03/
      data t_data_old( 44) /   0.613808411719D+03/
      data t_data_old( 45) /   0.626751944545D+03/
      data t_data_old( 46) /   0.640021153223D+03/
      data t_data_old( 47) /   0.646778108269D+03/
      data t_data_old( 48) /   0.653616475892D+03/
      data t_data_old( 49) /   0.660535909694D+03/
      data t_data_old( 50) /   0.667536042128D+03/
      data t_data_old( 51) /   0.674616484717D+03/
      data t_data_old( 52) /   0.681776828247D+03/
      data t_data_old( 53) /   0.689016642934D+03/
      data t_data_old( 54) /   0.696335478569D+03/
      data t_data_old( 55) /   0.703732864616D+03/
      data t_data_old( 56) /   0.711208310280D+03/
      data t_data_old( 57) /   0.718761304525D+03/
      data t_data_old( 58) /   0.726391316044D+03/
      data t_data_old( 59) /   0.730235124626D+03/
      data t_data_old( 60) /   0.734098074953D+03/
      data t_data_old( 61) /   0.737980094205D+03/
      data t_data_old( 62) /   0.741881108249D+03/
      data t_data_old( 63) /   0.749739642169D+03/
      data t_data_old( 64) /   0.757672863091D+03/
      data t_data_old( 65) /   0.765680135304D+03/
      data t_data_old( 66) /   0.773760801178D+03/
      data t_data_old( 67) /   0.781914180580D+03/
      data t_data_old( 68) /   0.790139570186D+03/
      data t_data_old( 69) /   0.798436242688D+03/
      data t_data_old( 70) /   0.806803445880D+03/
      data t_data_old( 71) /   0.815240401625D+03/
      data t_data_old( 72) /   0.823746304702D+03/
      data t_data_old( 73) /   0.832320321528D+03/
      data t_data_old( 74) /   0.840961588750D+03/
      data t_data_old( 75) /   0.849669211710D+03/
      data t_data_old( 76) /   0.858442262774D+03/
      data t_data_old( 77) /   0.876180129380D+03/
      data t_data_old( 78) /   0.894165517659D+03/
      data t_data_old( 79) /   0.912389760244D+03/
      data t_data_old( 80) /   0.930843478703D+03/
      data t_data_old( 81) /   0.949516507809D+03/
      data t_data_old( 82) /   0.968397815750D+03/
      data t_data_old( 83) /   0.987475421670D+03/
      data t_data_old( 84) /   0.100673631209D+04/
      data t_data_old( 85) /   0.102616628379D+04/
      data t_data_old( 86) /   0.103593929029D+04/
      data t_data_old( 87) /   0.104574892620D+04/
      data t_data_old( 88) /   0.106547088379D+04/
      data t_data_old( 89) /   0.107537667398D+04/
      data t_data_old( 90) /   0.108530999416D+04/
      data t_data_old( 91) /   0.109526841725D+04/
      data t_data_old( 92) /   0.110524932240D+04/
      data t_data_old( 93) /   0.111525004872D+04/
      data t_data_old( 94) /   0.112526784615D+04/
      data t_data_old( 95) /   0.113529987483D+04/
      data t_data_old( 96) /   0.114534320467D+04/
      data t_data_old( 97) /   0.115539481518D+04/
      data t_data_old( 98) /   0.116545159561D+04/
      data t_data_old( 99) /   0.117551034535D+04/
      data t_data_old(100) /   0.118556777459D+04/
      data t_data_old(101) /   0.119562050535D+04/
      data t_data_old(102) /   0.120566507280D+04/
      data t_data_old(103) /   0.121569792689D+04/
      data t_data_old(104) /   0.122571543436D+04/
      data t_data_old(105) /   0.123571388107D+04/
      data t_data_old(106) /   0.124568947478D+04/
      data t_data_old(107) /   0.125563834818D+04/
      data t_data_old(108) /   0.126555656244D+04/
      data t_data_old(109) /   0.127544011113D+04/
      data t_data_old(110) /   0.128528492453D+04/
      data t_data_old(111) /   0.129508687447D+04/
      data t_data_old(112) /   0.130484177953D+04/
      data t_data_old(113) /   0.131454541082D+04/
      data t_data_old(114) /   0.132419349812D+04/
      data t_data_old(115) /   0.133235052449D+04/
      data t_data_old(116) /   0.134046156686D+04/
      data t_data_old(117) /   0.134852394906D+04/
      data t_data_old(118) /   0.135653498754D+04/
      data t_data_old(119) /   0.136449199529D+04/
      data t_data_old(120) /   0.137239228577D+04/
      data t_data_old(121) /   0.138023317719D+04/
      data t_data_old(122) /   0.138801199688D+04/
      data t_data_old(123) /   0.139572608589D+04/
      data t_data_old(124) /   0.140337280379D+04/
      data t_data_old(125) /   0.141094953355D+04/
      data t_data_old(126) /   0.141845368673D+04/
      data t_data_old(127) /   0.142588270865D+04/
      data t_data_old(128) /   0.143323408384D+04/
      data t_data_old(129) /   0.144050534153D+04/
      data t_data_old(130) /   0.144769406128D+04/
      data t_data_old(131) /   0.145479787862D+04/
      data t_data_old(132) /   0.146181449088D+04/
      data t_data_old(133) /   0.146874166290D+04/
      data t_data_old(134) /   0.147557723289D+04/
      data t_data_old(135) /   0.148231911815D+04/
      data t_data_old(136) /   0.148896532085D+04/
      data t_data_old(137) /   0.149551393366D+04/
      data t_data_old(138) /   0.150196314527D+04/
      data t_data_old(139) /   0.150831124583D+04/
      data t_data_old(140) /   0.151455663215D+04/
      data t_data_old(141) /   0.152069781270D+04/
      data t_data_old(142) /   0.152673341238D+04/
      data t_data_old(143) /   0.153266217702D+04/
      data t_data_old(144) /   0.153848297754D+04/
      data t_data_old(145) /   0.154419481377D+04/
      data t_data_old(146) /   0.154979681787D+04/
      data t_data_old(147) /   0.155528825740D+04/
      data t_data_old(148) /   0.156066853786D+04/
      data t_data_old(149) /   0.156593720487D+04/
      data t_data_old(150) /   0.157109394570D+04/
      data t_data_old(151) /   0.157613859048D+04/
      data t_data_old(152) /   0.158107111269D+04/
      data t_data_old(153) /   0.158589162922D+04/
      data t_data_old(154) /   0.159060039984D+04/
      data t_data_old(155) /   0.159519782605D+04/
      data t_data_old(156) /   0.159968444944D+04/
      data t_data_old(157) /   0.160406094943D+04/
      data t_data_old(158) /   0.160832814046D+04/
      data t_data_old(159) /   0.161248696859D+04/
      data t_data_old(160) /   0.161653850766D+04/
      data t_data_old(161) /   0.162048395485D+04/
      data t_data_old(162) /   0.162432462576D+04/
      data t_data_old(163) /   0.162806194910D+04/
      data t_data_old(164) /   0.163169746083D+04/
      data t_data_old(165) /   0.163523279806D+04/
      data t_data_old(166) /   0.163866969250D+04/
      data t_data_old(167) /   0.164200996365D+04/
      data t_data_old(168) /   0.164525551170D+04/
      data t_data_old(169) /   0.164840831029D+04/
      data t_data_old(170) /   0.165147039903D+04/
      data t_data_old(171) /   0.165444387595D+04/
      data t_data_old(172) /   0.165733088991D+04/
      data t_data_old(173) /   0.166285247754D+04/
      data t_data_old(174) /   0.166805645889D+04/
      data t_data_old(175) /   0.167296105346D+04/
      data t_data_old(176) /   0.167758460573D+04/
      data t_data_old(177) /   0.168194537970D+04/
      data t_data_old(178) /   0.168606137252D+04/
      data t_data_old(179) /   0.168995015032D+04/
      data t_data_old(180) /   0.169362870800D+04/
      data t_data_old(181) /   0.169711335356D+04/
      data t_data_old(182) /   0.170041961704D+04/
      data t_data_old(183) /   0.170356218304D+04/
      data t_data_old(184) /   0.170655484528D+04/
      data t_data_old(185) /   0.170941048125D+04/
      data t_data_old(186) /   0.171214104456D+04/
      data t_data_old(187) /   0.171475757256D+04/
      data t_data_old(188) /   0.171727020667D+04/
      data t_data_old(189) /   0.171968822300D+04/
      data t_data_old(190) /   0.172202007088D+04/
      data t_data_old(191) /   0.172427341725D+04/
      data t_data_old(192) /   0.172645519488D+04/
      data t_data_old(193) /   0.172857165300D+04/
      data t_data_old(194) /   0.173062840876D+04/
      data t_data_old(195) /   0.173263049847D+04/
      data t_data_old(196) /   0.173458242772D+04/
      data t_data_old(197) /   0.173648821967D+04/
      data t_data_old(198) /   0.173835146096D+04/
      data t_data_old(199) /   0.174017534504D+04/
      data t_data_old(200) /   0.174196271243D+04/
      data t_data_old(201) /   0.174371608810D+04/
      data t_data_old(202) /   0.174712768517D+04/
      data t_data_old(203) /   0.175042823386D+04/
      data t_data_old(204) /   0.175362944740D+04/
      data t_data_old(205) /   0.175674063688D+04/
      data t_data_old(206) /   0.175976925308D+04/
      data t_data_old(207) /   0.176272131286D+04/
      data t_data_old(208) /   0.176560173150D+04/
      data t_data_old(209) /   0.176841457948D+04/
      data t_data_old(210) /   0.177116327949D+04/
      data t_data_old(211) /   0.177385075692D+04/
      data t_data_old(212) /   0.177647955425D+04/
      data t_data_old(213) /   0.178156807902D+04/
      data t_data_old(214) /   0.178644987684D+04/
      data t_data_old(215) /   0.179113775639D+04/
      data t_data_old(216) /   0.179564296780D+04/
      data t_data_old(217) /   0.180413685223D+04/
      data t_data_old(218) /   0.181203310087D+04/
      data t_data_old(219) /   0.181939231931D+04/
      data t_data_old(220) /   0.182626722557D+04/
      data t_data_old(221) /   0.183270383545D+04/
      data t_data_old(222) /   0.183874239939D+04/
      data t_data_old(223) /   0.184441819200D+04/
      data t_data_old(224) /   0.184976218718D+04/
      data t_data_old(225) /   0.185480163413D+04/
      data t_data_old(226) /   0.185956054498D+04/
      data t_data_old(227) /   0.186406010256D+04/
      data t_data_old(228) /   0.186831898938D+04/
      data t_data_old(229) /   0.187616494362D+04/
      data t_data_old(230) /   0.188327868997D+04/
      data t_data_old(231) /   0.188975358800D+04/
      data t_data_old(232) /   0.189566641052D+04/
      data t_data_old(233) /   0.190108051131D+04/
      data t_data_old(234) /   0.190604807068D+04/
      data t_data_old(235) /   0.191477176491D+04/
      data t_data_old(236) /   0.192230864634D+04/
      data t_data_old(237) /   0.192887122215D+04/
      data t_data_old(238) /   0.193462083295D+04/
      data t_data_old(239) /   0.193967931204D+04/
      data t_data_old(240) /   0.194798796845D+04/
      data t_data_old(241) /   0.195474979168D+04/
      data t_data_old(242) /   0.196033679375D+04/
      data t_data_old(243) /   0.196500901286D+04/
      data t_data_old(244) /   0.196894759839D+04/
      data t_data_old(245) /   0.197496994536D+04/
      data t_data_old(246) /   0.197959444603D+04/
      data t_data_old(247) /   0.198322498771D+04/
      data t_data_old(248) /   0.198611993684D+04/
      data t_data_old(249) /   0.199014258591D+04/
      data t_data_old(250) /   0.199278405680D+04/
      data t_data_old(251) /   0.199278405680D+04/

c      do i=1,N_old
c         print *,x_data_old(i),T_data_old(i)
c      enddo
c      stop

      debug = .true.
      debug = .false.
      !
      ! Given T coming in, find X corresponding to T.
      !
      if (debug) print*, 't = ', t
      if (t .le. t_data_old(1)) then
         x = x_data_old(1)
      else if (t .ge. t_data_old(N_old)) then
         x = x_data_old(N_old)
      else
         do i = 2, N_old
            if (t .le. t_data_old(i)) exit
         end do
         if (debug) then
          print*, 't: [', t_data_old(i-1), ' , ', t_data_old(i), ']'
          print*, 'x: [', x_data_old(i-1), ' , ', x_data_old(i), ']'
         end if
         frac =(t-t_data_old(i-1))/(t_data_old(i)-t_data_old(i-1))
         x = x_data_old(i-1)+frac*(x_data_old(i)-x_data_old(i-1))
      end if

      if (debug) print*, 'x = ', x
      !
      ! Shift x to deal with differences in peak heat release.
      !
      x = phr_new+(x-phr_old)

      if (debug) print*, 'x in new coordinates: ', x

      call pmf_C3H8_07(x, x, tmp, Nvals_new)
      if (Nvals_new.gt.Nvals_newMAX) then
         print *,'convert: tmp array too small for pmf vals'
         call bl_abort()
      endif      
      Nspec_new = Nvals_new - 3
      !
      ! HACK
      !
      !if (tmp(1) .gt. 1850.0D0) then
      !   call pmf_hflame27(x_hi_new, x_hi_new, tmp)
      !endif
      !
      ! Store temperature after mass fractions.
      !
      y(Nspec_new+1) = tmp(1)
      if (debug) print *,'new temperature: ',tmp(1)

      do i = 1, Nspec_new
         y(i) = tmp(3+i)
      end do

      sum = 0.0d0
      do i = 1, Nspec_new
         y(i) = MAX(y(i), 0.0d0)
         sum = sum + y(i)
      enddo
      y(Nspec_new) = y(Nspec_new) + 1.0d0 - sum

      ! Note that this function returns mole fractions
      end


      integer function FORT_TfromHYpt_new(T,Hin,Y,errMax,NiterMAX,res)
      implicit none
      REAL_T T,Y(*),H,Hin
      REAL_T TMIN,TMAX,errMAX
      integer NiterMAX,Niter,n,NiterDAMP
      parameter (TMIN=250, TMAX=5000)
      REAL_T  T0,cp,cv,dH,temp,RoverWbar,Wbar,RU,RUC,P1ATM
      REAL_T res(0:NiterMAX-1),dT, Htarg
      logical out_of_bounds, converged, soln_bad, stalled
      REAL_T h300,cp300,h6500,cp6500
      integer ihitlo,ihithi,j
      REAL_T RWRK
      integer IWRK,Nelt,Nspec,Nreac,Nfit,isioproc

      out_of_bounds(temp) = (temp.lt.TMIN-one) .or. (temp.gt.TMAX)

      call bl_pd_is_ioproc(isioproc)
      call CKINDX_new(IWRK,RWRK,Nelt,Nspec,Nreac,Nfit)

      NiterDAMP = NiterMAX
      if ((T.GE.TMIN).and.(T.LE.TMAX)) then
         T0 = T
      else
         T0 = half*(TMIN+TMAX)
         T = T0
      end if
      Niter = 0
      dH = zero
      soln_bad = .FALSE.
      Htarg = Hin * 1.d4
      ihitlo = 0
      ihithi = 0

      CALL CKHBMS_new(T,Y,IWRK,RWRK,H)
      dH = two*ABS(H - Htarg)/(one + ABS(H) + ABS(Htarg))
      res(Niter) = dH
      converged = dH.le.errMAX

      do while ((.not.converged) .and. (.not.soln_bad))

         CALL CKCPBS_new(T,Y,IWRK,RWRK,cp)
         dT = (Htarg - H)/cp
         if ((Niter.le.NiterDAMP).and.(T+dT.ge.TMAX)) then
            T = TMAX
            ihithi = 1
         else if ((Niter.le.NiterDAMP).and.(T+dT.le.TMIN)) then
            T = TMIN
            ihitlo = 1
         else
            T = T + dT
         end if
         soln_bad = out_of_bounds(T)
         if (soln_bad) then
            FORT_TfromHYpt_new = -1
            goto 100
         else
            CALL CKHBMS_new(T,Y,IWRK,RWRK,H)
            dH = two*ABS(H - Htarg)/(one + ABS(H) + ABS(Htarg))
            res(Niter) = dH
            Niter = Niter + 1
         end if
         if (Niter .ge. NiterMAX) then
            FORT_TfromHYpt_new = -2
            goto 100
         endif
         converged = (dH.le.errMAX) .or. (ABS(dT).le.errMAX)

         if ((ihitlo.eq.1).and.(H.gt.Htarg)) then
            T = TMIN
            CALL CKHBMS_new(T,Y,IWRK,RWRK,h300)
            CALL CKCPBS_new(T,Y,IWRK,RWRK,cp300)
            T=TMIN+(Htarg-h300)/cp300
            converged = .true.
         endif
         if ((ihithi.eq.1).and.(H.lt.Htarg)) then
            T = TMAX
            CALL CKHBMS_new(T,Y,IWRK,RWRK,h6500)
            CALL CKCPBS_new(T,Y,IWRK,RWRK,cp6500)
            T=TMAX+(Htarg-h6500)/cp6500
            converged = .true.
         endif
      end do
c
c     Set max iters taken during this solve and exit
c
      FORT_TfromHYpt_new = Niter
      return
c
c     Error condition....dump state and bail out
c
 100  continue

      if (isioproc.eq.1) then
      write(6,997) 'T from (H,Y) new: failed'
      write(6,997) 'iterations tried = ',Niter
      write(6,998) 'initial T = ',T0
      write(6,998) 'current T = ',T
      write(6,*) 'species mass fracs  1:',Nspec
      do n = 1,Nspec
         write(6,998) '  ',Y(n)
      end do
      write(6,998)
      write(6,998) 'residual:'
      do n = 0,Niter-1
         write(6,998) '  ',res(n)
      end do

 997  format(a,3(i4,a))
 998  format(a,d21.12)
      endif
      end

c      program main
c      implicit none
c      integer N,i
c      double precision x, Y(50), XLO, XHI
c      parameter (N=101)
c      XLO = 298.d0
c      XHI = 2005.d0
cc      XLO = 0.d0
cc      XHI = 1.5d0
c      do i=1,N
c         x = XLO + (XHI-XLO) * (i-1)/(N-1)
cc         call pmf_hflame27(x,x,Y)
cc         print *,x,Y(1)
c         call FORT_CONVERT_DRM19_TO_H2 (x, Y)
c         print *,x,Y(10)
c      enddo
cc      call FORT_CONVERT_DRM19_TO_H2 (x, Y)
c      end       

