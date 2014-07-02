namespace LumiFilter
{

bool pass(unsigned int run, unsigned int lumi)
{
  // Fallback for MC and 2011 data
  if (run < 190000) return true;

  // 22Jan Re-Reco v1 JSON
  if (run == 198050 && (lumi >= 2 && lumi <= 155))
    return true;
  if (run == 201729 && (lumi >= 6 && lumi <= 20))
    return true;
  if (run == 201729 && (lumi >= 22 && lumi <= 75))
    return true;
  if (run == 201729 && (lumi >= 77 && lumi <= 126))
    return true;
  if (run == 201729 && (lumi >= 129 && lumi <= 154))
    return true;
  if (run == 201729 && (lumi >= 156 && lumi <= 216))
    return true;
  if (run == 201729 && (lumi >= 219 && lumi <= 244))
    return true;
  if (run == 201727 && (lumi >= 67 && lumi <= 185))
    return true;
  if (run == 201802 && (lumi >= 68 && lumi <= 209))
    return true;
  if (run == 201802 && (lumi >= 211 && lumi <= 214))
    return true;
  if (run == 201802 && (lumi >= 216 && lumi <= 220))
    return true;
  if (run == 201802 && (lumi >= 223 && lumi <= 288))
    return true;
  if (run == 201802 && (lumi >= 290 && lumi <= 296))
    return true;
  if (run == 199573 && (lumi >= 1 && lumi <= 22))
    return true;
  if (run == 199572 && (lumi >= 1 && lumi <= 317))
    return true;
  if (run == 199571 && (lumi >= 1 && lumi <= 184))
    return true;
  if (run == 199571 && (lumi >= 186 && lumi <= 360))
    return true;
  if (run == 199571 && (lumi >= 363 && lumi <= 561))
    return true;
  if (run == 199570 && (lumi >= 1 && lumi <= 17))
    return true;
  if (run == 199574 && (lumi >= 1 && lumi <= 53))
    return true;
  if (run == 199574 && (lumi >= 56 && lumi <= 153))
    return true;
  if (run == 199574 && (lumi >= 156 && lumi <= 246))
    return true;
  if (run == 195378 && (lumi >= 90 && lumi <= 117))
    return true;
  if (run == 195378 && (lumi >= 120 && lumi <= 127))
    return true;
  if (run == 195378 && (lumi >= 130 && lumi <= 185))
    return true;
  if (run == 195378 && (lumi >= 187 && lumi <= 204))
    return true;
  if (run == 195378 && (lumi >= 206 && lumi <= 302))
    return true;
  if (run == 195378 && (lumi >= 305 && lumi <= 542))
    return true;
  if (run == 195378 && (lumi >= 544 && lumi <= 565))
    return true;
  if (run == 195378 && (lumi >= 567 && lumi <= 645))
    return true;
  if (run == 195378 && (lumi >= 647 && lumi <= 701))
    return true;
  if (run == 195378 && (lumi >= 703 && lumi <= 734))
    return true;
  if (run == 195378 && (lumi >= 737 && lumi <= 1120))
    return true;
  if (run == 195378 && (lumi >= 1122 && lumi <= 1133))
    return true;
  if (run == 206401 && (lumi >= 69 && lumi <= 90))
    return true;
  if (run == 206401 && (lumi >= 92 && lumi <= 194))
    return true;
  if (run == 206401 && (lumi >= 197 && lumi <= 210))
    return true;
  if (run == 206401 && (lumi >= 212 && lumi <= 249))
    return true;
  if (run == 206401 && (lumi >= 251 && lumi <= 265))
    return true;
  if (run == 206401 && (lumi >= 267 && lumi <= 409))
    return true;
  if (run == 193093 && (lumi >= 1 && lumi <= 33))
    return true;
  if (run == 204577 && (lumi >= 1 && lumi <= 46))
    return true;
  if (run == 204577 && (lumi >= 49 && lumi <= 64))
    return true;
  if (run == 204577 && (lumi >= 67 && lumi <= 105))
    return true;
  if (run == 204577 && (lumi >= 107 && lumi <= 170))
    return true;
  if (run == 204577 && (lumi >= 173 && lumi <= 181))
    return true;
  if (run == 204577 && (lumi >= 183 && lumi <= 193))
    return true;
  if (run == 204577 && (lumi >= 196 && lumi <= 653))
    return true;
  if (run == 204577 && (lumi >= 656 && lumi <= 669))
    return true;
  if (run == 204577 && (lumi >= 671 && lumi <= 740))
    return true;
  if (run == 204577 && (lumi >= 742 && lumi <= 913))
    return true;
  if (run == 204577 && (lumi >= 915 && lumi <= 1057))
    return true;
  if (run == 204577 && (lumi >= 1059 && lumi <= 1115))
    return true;
  if (run == 204577 && (lumi >= 1117 && lumi <= 1282))
    return true;
  if (run == 205111 && (lumi >= 88 && lumi <= 390))
    return true;
  if (run == 205111 && (lumi >= 392 && lumi <= 441))
    return true;
  if (run == 205111 && (lumi >= 444 && lumi <= 446))
    return true;
  if (run == 198954 && (lumi >= 108 && lumi <= 156))
    return true;
  if (run == 198954 && (lumi >= 159 && lumi <= 277))
    return true;
  if (run == 198955 && (lumi >= 1 && lumi <= 45))
    return true;
  if (run == 198955 && (lumi >= 47 && lumi <= 50))
    return true;
  if (run == 198955 && (lumi >= 53 && lumi <= 220))
    return true;
  if (run == 198955 && (lumi >= 223 && lumi <= 269))
    return true;
  if (run == 198955 && (lumi >= 271 && lumi <= 284))
    return true;
  if (run == 198955 && (lumi >= 286 && lumi <= 338))
    return true;
  if (run == 198955 && (lumi >= 340 && lumi <= 580))
    return true;
  if (run == 198955 && (lumi >= 583 && lumi <= 742))
    return true;
  if (run == 198955 && (lumi >= 744 && lumi <= 910))
    return true;
  if (run == 198955 && (lumi >= 913 && lumi <= 946))
    return true;
  if (run == 198955 && (lumi >= 949 && lumi <= 1162))
    return true;
  if (run == 198955 && (lumi >= 1165 && lumi <= 1169))
    return true;
  if (run == 198955 && (lumi >= 1172 && lumi <= 1182))
    return true;
  if (run == 198955 && (lumi >= 1185 && lumi <= 1188))
    return true;
  if (run == 198955 && (lumi >= 1190 && lumi <= 1246))
    return true;
  if (run == 198955 && (lumi >= 1249 && lumi <= 1304))
    return true;
  if (run == 198955 && (lumi >= 1306 && lumi <= 1467))
    return true;
  if (run == 198955 && (lumi >= 1470 && lumi <= 1485))
    return true;
  if (run == 198955 && (lumi >= 1487 && lumi <= 1552))
    return true;
  if (run == 195109 && (lumi >= 98 && lumi <= 241))
    return true;
  if (run == 202060 && (lumi >= 76 && lumi <= 142))
    return true;
  if (run == 202060 && (lumi >= 144 && lumi <= 154))
    return true;
  if (run == 202060 && (lumi >= 156 && lumi <= 244))
    return true;
  if (run == 202060 && (lumi >= 246 && lumi <= 497))
    return true;
  if (run == 202060 && (lumi >= 499 && lumi <= 642))
    return true;
  if (run == 202060 && (lumi >= 644 && lumi <= 682))
    return true;
  if (run == 202060 && (lumi >= 684 && lumi <= 743))
    return true;
  if (run == 202060 && (lumi >= 746 && lumi <= 936))
    return true;
  if (run == 194108 && (lumi >= 81 && lumi <= 161))
    return true;
  if (run == 194108 && (lumi >= 164 && lumi <= 264))
    return true;
  if (run == 194108 && (lumi >= 266 && lumi <= 373))
    return true;
  if (run == 194108 && (lumi >= 376 && lumi <= 396))
    return true;
  if (run == 194108 && (lumi >= 398 && lumi <= 433))
    return true;
  if (run == 194108 && (lumi >= 436 && lumi <= 452))
    return true;
  if (run == 194108 && (lumi >= 454 && lumi <= 577))
    return true;
  if (run == 194108 && (lumi >= 579 && lumi <= 590))
    return true;
  if (run == 194108 && (lumi >= 593 && lumi <= 668))
    return true;
  if (run == 194108 && (lumi >= 671 && lumi <= 872))
    return true;
  if (run == 200244 && (lumi >= 3 && lumi <= 304))
    return true;
  if (run == 200244 && (lumi >= 307 && lumi <= 442))
    return true;
  if (run == 200244 && (lumi >= 445 && lumi <= 507))
    return true;
  if (run == 200244 && (lumi >= 510 && lumi <= 619))
    return true;
  if (run == 200245 && (lumi >= 1 && lumi <= 103))
    return true;
  if (run == 200245 && (lumi >= 105 && lumi <= 128))
    return true;
  if (run == 200245 && (lumi >= 131 && lumi <= 248))
    return true;
  if (run == 200245 && (lumi >= 251 && lumi <= 357))
    return true;
  if (run == 200243 && (lumi >= 69 && lumi <= 103))
    return true;
  if (run == 200243 && (lumi >= 106 && lumi <= 139))
    return true;
  if (run == 194778 && (lumi >= 60 && lumi <= 118))
    return true;
  if (run == 194778 && (lumi >= 120 && lumi <= 219))
    return true;
  if (run == 206331 && (lumi >= 91 && lumi <= 222))
    return true;
  if (run == 206331 && (lumi >= 225 && lumi <= 312))
    return true;
  if (run == 205236 && (lumi >= 1 && lumi <= 190))
    return true;
  if (run == 205236 && (lumi >= 193 && lumi <= 207))
    return true;
  if (run == 205236 && (lumi >= 209 && lumi <= 260))
    return true;
  if (run == 205236 && (lumi >= 263 && lumi <= 331))
    return true;
  if (run == 205236 && (lumi >= 334 && lumi <= 352))
    return true;
  if (run == 194479 && (lumi >= 1 && lumi <= 44))
    return true;
  if (run == 194479 && (lumi >= 165 && lumi <= 232))
    return true;
  if (run == 194479 && (lumi >= 235 && lumi <= 262))
    return true;
  if (run == 194479 && (lumi >= 265 && lumi <= 374))
    return true;
  if (run == 194479 && (lumi >= 377 && lumi <= 431))
    return true;
  if (run == 194479 && (lumi >= 434 && lumi <= 489))
    return true;
  if (run == 194479 && (lumi >= 492 && lumi <= 529))
    return true;
  if (run == 194479 && (lumi >= 531 && lumi <= 566))
    return true;
  if (run == 205233 && (lumi >= 94 && lumi <= 153))
    return true;
  if (run == 203912 && (lumi >= 1 && lumi <= 306))
    return true;
  if (run == 203912 && (lumi >= 308 && lumi <= 566))
    return true;
  if (run == 203912 && (lumi >= 569 && lumi <= 609))
    return true;
  if (run == 203912 && (lumi >= 611 && lumi <= 698))
    return true;
  if (run == 203912 && (lumi >= 701 && lumi <= 820))
    return true;
  if (run == 203912 && (lumi >= 823 && lumi <= 865))
    return true;
  if (run == 203912 && (lumi >= 867 && lumi <= 1033))
    return true;
  if (run == 203912 && (lumi >= 1035 && lumi <= 1321))
    return true;
  if (run == 206102 && (lumi >= 83 && lumi <= 116))
    return true;
  if (run == 206102 && (lumi >= 120 && lumi <= 130))
    return true;
  if (run == 206102 && (lumi >= 133 && lumi <= 208))
    return true;
  if (run == 206102 && (lumi >= 211 && lumi <= 235))
    return true;
  if (run == 206102 && (lumi >= 238 && lumi <= 246))
    return true;
  if (run == 206102 && (lumi >= 249 && lumi <= 278))
    return true;
  if (run == 206102 && (lumi >= 281 && lumi <= 349))
    return true;
  if (run == 205238 && (lumi >= 1 && lumi <= 6))
    return true;
  if (run == 205238 && (lumi >= 9 && lumi <= 199))
    return true;
  if (run == 205238 && (lumi >= 202 && lumi <= 254))
    return true;
  if (run == 205238 && (lumi >= 256 && lumi <= 304))
    return true;
  if (run == 205238 && (lumi >= 306 && lumi <= 355))
    return true;
  if (run == 205238 && (lumi >= 358 && lumi <= 381))
    return true;
  if (run == 205238 && (lumi >= 384 && lumi <= 596))
    return true;
  if (run == 205238 && (lumi >= 598 && lumi <= 617))
    return true;
  if (run == 193575 && (lumi >= 48 && lumi <= 173))
    return true;
  if (run == 193575 && (lumi >= 176 && lumi <= 349))
    return true;
  if (run == 193575 && (lumi >= 351 && lumi <= 394))
    return true;
  if (run == 193575 && (lumi >= 397 && lumi <= 415))
    return true;
  if (run == 193575 && (lumi >= 417 && lumi <= 658))
    return true;
  if (run == 193575 && (lumi >= 660 && lumi <= 752))
    return true;
  if (run == 191264 && (lumi >= 59 && lumi <= 79))
    return true;
  if (run == 191264 && (lumi >= 82 && lumi <= 152))
    return true;
  if (run == 191264 && (lumi >= 155 && lumi <= 189))
    return true;
  if (run == 201229 && (lumi >= 1 && lumi <= 5))
    return true;
  if (run == 201229 && (lumi >= 8 && lumi <= 26))
    return true;
  if (run == 201229 && (lumi >= 29 && lumi <= 73))
    return true;
  if (run == 199739 && (lumi >= 66 && lumi <= 133))
    return true;
  if (run == 195950 && (lumi >= 1 && lumi <= 71))
    return true;
  if (run == 195950 && (lumi >= 73 && lumi <= 138))
    return true;
  if (run == 195950 && (lumi >= 141 && lumi <= 169))
    return true;
  if (run == 195950 && (lumi >= 172 && lumi <= 332))
    return true;
  if (run == 195950 && (lumi >= 335 && lumi <= 350))
    return true;
  if (run == 195950 && (lumi >= 353 && lumi <= 382))
    return true;
  if (run == 195950 && (lumi >= 385 && lumi <= 421))
    return true;
  if (run == 195950 && (lumi >= 424 && lumi <= 450))
    return true;
  if (run == 195950 && (lumi >= 453 && lumi <= 483))
    return true;
  if (run == 195950 && (lumi >= 485 && lumi <= 616))
    return true;
  if (run == 195950 && (lumi >= 619 && lumi <= 715))
    return true;
  if (run == 195950 && (lumi >= 718 && lumi <= 787))
    return true;
  if (run == 195950 && (lumi >= 789 && lumi <= 800))
    return true;
  if (run == 195950 && (lumi >= 803 && lumi <= 829))
    return true;
  if (run == 195950 && (lumi >= 831 && lumi <= 831))
    return true;
  if (run == 195950 && (lumi >= 833 && lumi <= 1587))
    return true;
  if (run == 201824 && (lumi >= 1 && lumi <= 139))
    return true;
  if (run == 201824 && (lumi >= 141 && lumi <= 176))
    return true;
  if (run == 201824 && (lumi >= 179 && lumi <= 286))
    return true;
  if (run == 201824 && (lumi >= 289 && lumi <= 492))
    return true;
  if (run == 199804 && (lumi >= 78 && lumi <= 88))
    return true;
  if (run == 199804 && (lumi >= 90 && lumi <= 181))
    return true;
  if (run == 199804 && (lumi >= 183 && lumi <= 235))
    return true;
  if (run == 199804 && (lumi >= 238 && lumi <= 278))
    return true;
  if (run == 199804 && (lumi >= 281 && lumi <= 290))
    return true;
  if (run == 199804 && (lumi >= 292 && lumi <= 519))
    return true;
  if (run == 199804 && (lumi >= 522 && lumi <= 575))
    return true;
  if (run == 199804 && (lumi >= 577 && lumi <= 628))
    return true;
  if (run == 199804 && (lumi >= 631 && lumi <= 632))
    return true;
  if (run == 191046 && (lumi >= 1 && lumi <= 21))
    return true;
  if (run == 191046 && (lumi >= 24 && lumi <= 82))
    return true;
  if (run == 191046 && (lumi >= 84 && lumi <= 88))
    return true;
  if (run == 191046 && (lumi >= 92 && lumi <= 116))
    return true;
  if (run == 191046 && (lumi >= 119 && lumi <= 180))
    return true;
  if (run == 191046 && (lumi >= 183 && lumi <= 183))
    return true;
  if (run == 191046 && (lumi >= 185 && lumi <= 239))
    return true;
  if (run == 191043 && (lumi >= 45 && lumi <= 46))
    return true;
  if (run == 206208 && (lumi >= 1 && lumi <= 470))
    return true;
  if (run == 206208 && (lumi >= 472 && lumi <= 518))
    return true;
  if (run == 207883 && (lumi >= 1 && lumi <= 1))
    return true;
  if (run == 207883 && (lumi >= 3 && lumi <= 4))
    return true;
  if (run == 207883 && (lumi >= 7 && lumi <= 75))
    return true;
  if (run == 196359 && (lumi >= 1 && lumi <= 2))
    return true;
  if (run == 207885 && (lumi >= 1 && lumi <= 90))
    return true;
  if (run == 207884 && (lumi >= 1 && lumi <= 106))
    return true;
  if (run == 207884 && (lumi >= 108 && lumi <= 183))
    return true;
  if (run == 200091 && (lumi >= 67 && lumi <= 67))
    return true;
  if (run == 200091 && (lumi >= 70 && lumi <= 151))
    return true;
  if (run == 200091 && (lumi >= 154 && lumi <= 172))
    return true;
  if (run == 200091 && (lumi >= 174 && lumi <= 187))
    return true;
  if (run == 200091 && (lumi >= 190 && lumi <= 196))
    return true;
  if (run == 200091 && (lumi >= 199 && lumi <= 201))
    return true;
  if (run == 200091 && (lumi >= 204 && lumi <= 425))
    return true;
  if (run == 200091 && (lumi >= 428 && lumi <= 535))
    return true;
  if (run == 200091 && (lumi >= 537 && lumi <= 607))
    return true;
  if (run == 200091 && (lumi >= 610 && lumi <= 879))
    return true;
  if (run == 200091 && (lumi >= 881 && lumi <= 943))
    return true;
  if (run == 200091 && (lumi >= 946 && lumi <= 999))
    return true;
  if (run == 200091 && (lumi >= 1001 && lumi <= 1025))
    return true;
  if (run == 200091 && (lumi >= 1027 && lumi <= 1132))
    return true;
  if (run == 200091 && (lumi >= 1135 && lumi <= 1339))
    return true;
  if (run == 200091 && (lumi >= 1341 && lumi <= 1433))
    return true;
  if (run == 200091 && (lumi >= 1435 && lumi <= 1450))
    return true;
  if (run == 200091 && (lumi >= 1453 && lumi <= 1523))
    return true;
  if (run == 200091 && (lumi >= 1526 && lumi <= 1664))
    return true;
  if (run == 200091 && (lumi >= 1667 && lumi <= 1680))
    return true;
  if (run == 200091 && (lumi >= 1683 && lumi <= 1710))
    return true;
  if (run == 196452 && (lumi >= 82 && lumi <= 112))
    return true;
  if (run == 196452 && (lumi >= 114 && lumi <= 490))
    return true;
  if (run == 196452 && (lumi >= 493 && lumi <= 586))
    return true;
  if (run == 196452 && (lumi >= 589 && lumi <= 618))
    return true;
  if (run == 196452 && (lumi >= 622 && lumi <= 668))
    return true;
  if (run == 196452 && (lumi >= 671 && lumi <= 716))
    return true;
  if (run == 196452 && (lumi >= 718 && lumi <= 726))
    return true;
  if (run == 196452 && (lumi >= 728 && lumi <= 956))
    return true;
  if (run == 196452 && (lumi >= 958 && lumi <= 1004))
    return true;
  if (run == 196452 && (lumi >= 1007 && lumi <= 1091))
    return true;
  if (run == 207889 && (lumi >= 1 && lumi <= 43))
    return true;
  if (run == 207889 && (lumi >= 47 && lumi <= 57))
    return true;
  if (run == 207889 && (lumi >= 60 && lumi <= 303))
    return true;
  if (run == 207889 && (lumi >= 306 && lumi <= 442))
    return true;
  if (run == 207889 && (lumi >= 445 && lumi <= 445))
    return true;
  if (run == 207889 && (lumi >= 447 && lumi <= 551))
    return true;
  if (run == 207889 && (lumi >= 553 && lumi <= 731))
    return true;
  if (run == 207889 && (lumi >= 733 && lumi <= 907))
    return true;
  if (run == 207889 && (lumi >= 910 && lumi <= 945))
    return true;
  if (run == 196357 && (lumi >= 1 && lumi <= 4))
    return true;
  if (run == 206207 && (lumi >= 82 && lumi <= 130))
    return true;
  if (run == 206207 && (lumi >= 132 && lumi <= 176))
    return true;
  if (run == 206207 && (lumi >= 179 && lumi <= 194))
    return true;
  if (run == 206207 && (lumi >= 196 && lumi <= 388))
    return true;
  if (run == 206207 && (lumi >= 390 && lumi <= 419))
    return true;
  if (run == 206207 && (lumi >= 422 && lumi <= 447))
    return true;
  if (run == 206207 && (lumi >= 450 && lumi <= 569))
    return true;
  if (run == 206207 && (lumi >= 572 && lumi <= 690))
    return true;
  if (run == 202973 && (lumi >= 1 && lumi <= 234))
    return true;
  if (run == 202973 && (lumi >= 237 && lumi <= 305))
    return true;
  if (run == 202973 && (lumi >= 308 && lumi <= 437))
    return true;
  if (run == 202973 && (lumi >= 439 && lumi <= 530))
    return true;
  if (run == 202973 && (lumi >= 532 && lumi <= 541))
    return true;
  if (run == 202973 && (lumi >= 544 && lumi <= 552))
    return true;
  if (run == 202973 && (lumi >= 555 && lumi <= 851))
    return true;
  if (run == 202973 && (lumi >= 853 && lumi <= 1408))
    return true;
  if (run == 200525 && (lumi >= 77 && lumi <= 149))
    return true;
  if (run == 200525 && (lumi >= 151 && lumi <= 164))
    return true;
  if (run == 200525 && (lumi >= 166 && lumi <= 190))
    return true;
  if (run == 200525 && (lumi >= 193 && lumi <= 276))
    return true;
  if (run == 200525 && (lumi >= 278 && lumi <= 311))
    return true;
  if (run == 200525 && (lumi >= 314 && lumi <= 464))
    return true;
  if (run == 200525 && (lumi >= 467 && lumi <= 488))
    return true;
  if (run == 200525 && (lumi >= 491 && lumi <= 674))
    return true;
  if (run == 200525 && (lumi >= 676 && lumi <= 704))
    return true;
  if (run == 200525 && (lumi >= 707 && lumi <= 755))
    return true;
  if (run == 200525 && (lumi >= 757 && lumi <= 895))
    return true;
  if (run == 200525 && (lumi >= 898 && lumi <= 937))
    return true;
  if (run == 200525 && (lumi >= 939 && lumi <= 990))
    return true;
  if (run == 205303 && (lumi >= 35 && lumi <= 54))
    return true;
  if (run == 205303 && (lumi >= 90 && lumi <= 132))
    return true;
  if (run == 205303 && (lumi >= 135 && lumi <= 144))
    return true;
  if (run == 199608 && (lumi >= 60 && lumi <= 157))
    return true;
  if (run == 199608 && (lumi >= 159 && lumi <= 209))
    return true;
  if (run == 199608 && (lumi >= 211 && lumi <= 341))
    return true;
  if (run == 199608 && (lumi >= 344 && lumi <= 390))
    return true;
  if (run == 199608 && (lumi >= 392 && lumi <= 461))
    return true;
  if (run == 199608 && (lumi >= 464 && lumi <= 800))
    return true;
  if (run == 199608 && (lumi >= 802 && lumi <= 1064))
    return true;
  if (run == 199608 && (lumi >= 1067 && lumi <= 1392))
    return true;
  if (run == 199608 && (lumi >= 1395 && lumi <= 1630))
    return true;
  if (run == 199608 && (lumi >= 1633 && lumi <= 1904))
    return true;
  if (run == 199608 && (lumi >= 1907 && lumi <= 1962))
    return true;
  if (run == 199608 && (lumi >= 1965 && lumi <= 2252))
    return true;
  if (run == 199608 && (lumi >= 2255 && lumi <= 2422))
    return true;
  if (run == 200188 && (lumi >= 1 && lumi <= 24))
    return true;
  if (run == 200188 && (lumi >= 27 && lumi <= 28))
    return true;
  if (run == 200188 && (lumi >= 31 && lumi <= 76))
    return true;
  if (run == 200188 && (lumi >= 79 && lumi <= 271))
    return true;
  if (run == 200188 && (lumi >= 274 && lumi <= 352))
    return true;
  if (run == 207233 && (lumi >= 1 && lumi <= 119))
    return true;
  if (run == 207233 && (lumi >= 121 && lumi <= 148))
    return true;
  if (run == 207231 && (lumi >= 70 && lumi <= 84))
    return true;
  if (run == 207231 && (lumi >= 86 && lumi <= 121))
    return true;
  if (run == 207231 && (lumi >= 123 && lumi <= 184))
    return true;
  if (run == 207231 && (lumi >= 187 && lumi <= 189))
    return true;
  if (run == 207231 && (lumi >= 192 && lumi <= 303))
    return true;
  if (run == 207231 && (lumi >= 306 && lumi <= 354))
    return true;
  if (run == 207231 && (lumi >= 357 && lumi <= 481))
    return true;
  if (run == 207231 && (lumi >= 484 && lumi <= 504))
    return true;
  if (run == 207231 && (lumi >= 508 && lumi <= 549))
    return true;
  if (run == 207231 && (lumi >= 552 && lumi <= 626))
    return true;
  if (run == 207231 && (lumi >= 628 && lumi <= 690))
    return true;
  if (run == 207231 && (lumi >= 693 && lumi <= 875))
    return true;
  if (run == 207231 && (lumi >= 878 && lumi <= 1000))
    return true;
  if (run == 207231 && (lumi >= 1003 && lumi <= 1170))
    return true;
  if (run == 207231 && (lumi >= 1173 && lumi <= 1187))
    return true;
  if (run == 207231 && (lumi >= 1189 && lumi <= 1227))
    return true;
  if (run == 207231 && (lumi >= 1229 && lumi <= 1415))
    return true;
  if (run == 207231 && (lumi >= 1418 && lumi <= 1445))
    return true;
  if (run == 207231 && (lumi >= 1447 && lumi <= 1505))
    return true;
  if (run == 206868 && (lumi >= 1 && lumi <= 3))
    return true;
  if (run == 206868 && (lumi >= 10 && lumi <= 16))
    return true;
  if (run == 206594 && (lumi >= 72 && lumi <= 107))
    return true;
  if (run == 206594 && (lumi >= 110 && lumi <= 246))
    return true;
  if (run == 206594 && (lumi >= 249 && lumi <= 281))
    return true;
  if (run == 206595 && (lumi >= 1 && lumi <= 34))
    return true;
  if (run == 206595 && (lumi >= 37 && lumi <= 42))
    return true;
  if (run == 206595 && (lumi >= 45 && lumi <= 193))
    return true;
  if (run == 206596 && (lumi >= 1 && lumi <= 13))
    return true;
  if (run == 206596 && (lumi >= 15 && lumi <= 220))
    return true;
  if (run == 206596 && (lumi >= 222 && lumi <= 228))
    return true;
  if (run == 206596 && (lumi >= 231 && lumi <= 236))
    return true;
  if (run == 206596 && (lumi >= 239 && lumi <= 292))
    return true;
  if (run == 206596 && (lumi >= 295 && lumi <= 695))
    return true;
  if (run == 206596 && (lumi >= 697 && lumi <= 728))
    return true;
  if (run == 206596 && (lumi >= 730 && lumi <= 810))
    return true;
  if (run == 206869 && (lumi >= 1 && lumi <= 251))
    return true;
  if (run == 206869 && (lumi >= 253 && lumi <= 271))
    return true;
  if (run == 206869 && (lumi >= 274 && lumi <= 502))
    return true;
  if (run == 206869 && (lumi >= 507 && lumi <= 520))
    return true;
  if (run == 206869 && (lumi >= 522 && lumi <= 566))
    return true;
  if (run == 206869 && (lumi >= 568 && lumi <= 752))
    return true;
  if (run == 206598 && (lumi >= 1 && lumi <= 81))
    return true;
  if (run == 206598 && (lumi >= 83 && lumi <= 103))
    return true;
  if (run == 206598 && (lumi >= 105 && lumi <= 588))
    return true;
  if (run == 206598 && (lumi >= 591 && lumi <= 657))
    return true;
  if (run == 206598 && (lumi >= 659 && lumi <= 719))
    return true;
  if (run == 199429 && (lumi >= 1 && lumi <= 28))
    return true;
  if (run == 199429 && (lumi >= 30 && lumi <= 36))
    return true;
  if (run == 199429 && (lumi >= 39 && lumi <= 55))
    return true;
  if (run == 199429 && (lumi >= 58 && lumi <= 101))
    return true;
  if (run == 199429 && (lumi >= 103 && lumi <= 148))
    return true;
  if (run == 199429 && (lumi >= 151 && lumi <= 154))
    return true;
  if (run == 199428 && (lumi >= 61 && lumi <= 197))
    return true;
  if (run == 199428 && (lumi >= 200 && lumi <= 210))
    return true;
  if (run == 199428 && (lumi >= 212 && lumi <= 382))
    return true;
  if (run == 199428 && (lumi >= 387 && lumi <= 414))
    return true;
  if (run == 199428 && (lumi >= 417 && lumi <= 436))
    return true;
  if (run == 199428 && (lumi >= 439 && lumi <= 530))
    return true;
  if (run == 199428 && (lumi >= 533 && lumi <= 648))
    return true;
  if (run == 200466 && (lumi >= 134 && lumi <= 274))
    return true;
  if (run == 196239 && (lumi >= 1 && lumi <= 59))
    return true;
  if (run == 196239 && (lumi >= 62 && lumi <= 154))
    return true;
  if (run == 196239 && (lumi >= 157 && lumi <= 272))
    return true;
  if (run == 196239 && (lumi >= 274 && lumi <= 373))
    return true;
  if (run == 196239 && (lumi >= 375 && lumi <= 432))
    return true;
  if (run == 196239 && (lumi >= 435 && lumi <= 465))
    return true;
  if (run == 196239 && (lumi >= 468 && lumi <= 647))
    return true;
  if (run == 196239 && (lumi >= 650 && lumi <= 706))
    return true;
  if (run == 196239 && (lumi >= 709 && lumi <= 1025))
    return true;
  if (run == 199876 && (lumi >= 1 && lumi <= 19))
    return true;
  if (run == 199876 && (lumi >= 22 && lumi <= 95))
    return true;
  if (run == 199876 && (lumi >= 97 && lumi <= 249))
    return true;
  if (run == 199876 && (lumi >= 252 && lumi <= 272))
    return true;
  if (run == 199876 && (lumi >= 274 && lumi <= 340))
    return true;
  if (run == 199876 && (lumi >= 343 && lumi <= 362))
    return true;
  if (run == 199876 && (lumi >= 365 && lumi <= 376))
    return true;
  if (run == 199877 && (lumi >= 1 && lumi <= 173))
    return true;
  if (run == 199877 && (lumi >= 175 && lumi <= 605))
    return true;
  if (run == 199877 && (lumi >= 607 && lumi <= 701))
    return true;
  if (run == 199877 && (lumi >= 703 && lumi <= 871))
    return true;
  if (run == 199875 && (lumi >= 70 && lumi <= 150))
    return true;
  if (run == 199875 && (lumi >= 152 && lumi <= 334))
    return true;
  if (run == 206446 && (lumi >= 92 && lumi <= 141))
    return true;
  if (run == 206446 && (lumi >= 143 && lumi <= 159))
    return true;
  if (run == 206446 && (lumi >= 162 && lumi <= 205))
    return true;
  if (run == 206446 && (lumi >= 208 && lumi <= 301))
    return true;
  if (run == 206446 && (lumi >= 304 && lumi <= 442))
    return true;
  if (run == 206446 && (lumi >= 445 && lumi <= 445))
    return true;
  if (run == 206446 && (lumi >= 448 && lumi <= 474))
    return true;
  if (run == 206446 && (lumi >= 476 && lumi <= 616))
    return true;
  if (run == 206446 && (lumi >= 619 && lumi <= 872))
    return true;
  if (run == 206446 && (lumi >= 874 && lumi <= 910))
    return true;
  if (run == 206446 && (lumi >= 912 && lumi <= 948))
    return true;
  if (run == 206446 && (lumi >= 950 && lumi <= 989))
    return true;
  if (run == 206446 && (lumi >= 992 && lumi <= 1030))
    return true;
  if (run == 206446 && (lumi >= 1033 && lumi <= 1075))
    return true;
  if (run == 206446 && (lumi >= 1109 && lumi <= 1149))
    return true;
  if (run == 190705 && (lumi >= 1 && lumi <= 5))
    return true;
  if (run == 190705 && (lumi >= 7 && lumi <= 65))
    return true;
  if (run == 190705 && (lumi >= 81 && lumi <= 336))
    return true;
  if (run == 190705 && (lumi >= 338 && lumi <= 350))
    return true;
  if (run == 190705 && (lumi >= 353 && lumi <= 383))
    return true;
  if (run == 206448 && (lumi >= 1 && lumi <= 143))
    return true;
  if (run == 206448 && (lumi >= 145 && lumi <= 559))
    return true;
  if (run == 206448 && (lumi >= 561 && lumi <= 1170))
    return true;
  if (run == 206448 && (lumi >= 1173 && lumi <= 1231))
    return true;
  if (run == 206448 && (lumi >= 1235 && lumi <= 1237))
    return true;
  if (run == 208551 && (lumi >= 119 && lumi <= 193))
    return true;
  if (run == 208551 && (lumi >= 195 && lumi <= 212))
    return true;
  if (run == 208551 && (lumi >= 215 && lumi <= 300))
    return true;
  if (run == 208551 && (lumi >= 303 && lumi <= 354))
    return true;
  if (run == 208551 && (lumi >= 356 && lumi <= 554))
    return true;
  if (run == 208551 && (lumi >= 557 && lumi <= 580))
    return true;
  if (run == 199973 && (lumi >= 73 && lumi <= 89))
    return true;
  if (run == 204544 && (lumi >= 1 && lumi <= 11))
    return true;
  if (run == 204544 && (lumi >= 13 && lumi <= 93))
    return true;
  if (run == 204544 && (lumi >= 96 && lumi <= 195))
    return true;
  if (run == 204544 && (lumi >= 197 && lumi <= 224))
    return true;
  if (run == 204544 && (lumi >= 226 && lumi <= 334))
    return true;
  if (run == 204544 && (lumi >= 337 && lumi <= 426))
    return true;
  if (run == 194691 && (lumi >= 61 && lumi <= 104))
    return true;
  if (run == 194691 && (lumi >= 107 && lumi <= 155))
    return true;
  if (run == 194691 && (lumi >= 158 && lumi <= 251))
    return true;
  if (run == 194691 && (lumi >= 254 && lumi <= 268))
    return true;
  if (run == 194691 && (lumi >= 271 && lumi <= 272))
    return true;
  if (run == 194691 && (lumi >= 275 && lumi <= 289))
    return true;
  if (run == 194691 && (lumi >= 292 && lumi <= 313))
    return true;
  if (run == 194912 && (lumi >= 53 && lumi <= 70))
    return true;
  if (run == 194912 && (lumi >= 72 && lumi <= 96))
    return true;
  if (run == 194912 && (lumi >= 98 && lumi <= 444))
    return true;
  if (run == 194912 && (lumi >= 446 && lumi <= 450))
    return true;
  if (run == 194912 && (lumi >= 453 && lumi <= 467))
    return true;
  if (run == 194912 && (lumi >= 470 && lumi <= 561))
    return true;
  if (run == 194912 && (lumi >= 564 && lumi <= 660))
    return true;
  if (run == 194912 && (lumi >= 663 && lumi <= 813))
    return true;
  if (run == 194912 && (lumi >= 815 && lumi <= 840))
    return true;
  if (run == 194912 && (lumi >= 843 && lumi <= 864))
    return true;
  if (run == 194912 && (lumi >= 866 && lumi <= 1004))
    return true;
  if (run == 194912 && (lumi >= 1007 && lumi <= 1025))
    return true;
  if (run == 194912 && (lumi >= 1027 && lumi <= 1067))
    return true;
  if (run == 194912 && (lumi >= 1069 && lumi <= 1137))
    return true;
  if (run == 194912 && (lumi >= 1140 && lumi <= 1166))
    return true;
  if (run == 194912 && (lumi >= 1168 && lumi <= 1249))
    return true;
  if (run == 194912 && (lumi >= 1251 && lumi <= 1304))
    return true;
  if (run == 194912 && (lumi >= 1307 && lumi <= 1444))
    return true;
  if (run == 194912 && (lumi >= 1447 && lumi <= 1487))
    return true;
  if (run == 194912 && (lumi >= 1489 && lumi <= 1503))
    return true;
  if (run == 194912 && (lumi >= 1506 && lumi <= 1662))
    return true;
  if (run == 194699 && (lumi >= 1 && lumi <= 30))
    return true;
  if (run == 194699 && (lumi >= 32 && lumi <= 52))
    return true;
  if (run == 194699 && (lumi >= 55 && lumi <= 64))
    return true;
  if (run == 194699 && (lumi >= 67 && lumi <= 71))
    return true;
  if (run == 194699 && (lumi >= 73 && lumi <= 154))
    return true;
  if (run == 194699 && (lumi >= 157 && lumi <= 215))
    return true;
  if (run == 194699 && (lumi >= 218 && lumi <= 238))
    return true;
  if (run == 194699 && (lumi >= 241 && lumi <= 259))
    return true;
  if (run == 194914 && (lumi >= 1 && lumi <= 38))
    return true;
  if (run == 194915 && (lumi >= 1 && lumi <= 74))
    return true;
  if (run == 202972 && (lumi >= 1 && lumi <= 30))
    return true;
  if (run == 202972 && (lumi >= 33 && lumi <= 184))
    return true;
  if (run == 202972 && (lumi >= 186 && lumi <= 290))
    return true;
  if (run == 202972 && (lumi >= 292 && lumi <= 295))
    return true;
  if (run == 202972 && (lumi >= 298 && lumi <= 371))
    return true;
  if (run == 202972 && (lumi >= 374 && lumi <= 429))
    return true;
  if (run == 202972 && (lumi >= 431 && lumi <= 544))
    return true;
  if (run == 190708 && (lumi >= 1 && lumi <= 189))
    return true;
  if (run == 198063 && (lumi >= 1 && lumi <= 37))
    return true;
  if (run == 198063 && (lumi >= 40 && lumi <= 72))
    return true;
  if (run == 198063 && (lumi >= 74 && lumi <= 124))
    return true;
  if (run == 198063 && (lumi >= 127 && lumi <= 294))
    return true;
  if (run == 199961 && (lumi >= 1 && lumi <= 211))
    return true;
  if (run == 199961 && (lumi >= 213 && lumi <= 287))
    return true;
  if (run == 202314 && (lumi >= 67 && lumi <= 104))
    return true;
  if (run == 202314 && (lumi >= 107 && lumi <= 265))
    return true;
  if (run == 202314 && (lumi >= 268 && lumi <= 278))
    return true;
  if (run == 196046 && (lumi >= 12 && lumi <= 40))
    return true;
  if (run == 196047 && (lumi >= 1 && lumi <= 64))
    return true;
  if (run == 196047 && (lumi >= 70 && lumi <= 75))
    return true;
  if (run == 200961 && (lumi >= 1 && lumi <= 115))
    return true;
  if (run == 196048 && (lumi >= 1 && lumi <= 44))
    return true;
  if (run == 196048 && (lumi >= 46 && lumi <= 48))
    return true;
  if (run == 190659 && (lumi >= 33 && lumi <= 167))
    return true;
  if (run == 201168 && (lumi >= 1 && lumi <= 37))
    return true;
  if (run == 201168 && (lumi >= 39 && lumi <= 275))
    return true;
  if (run == 201168 && (lumi >= 278 && lumi <= 481))
    return true;
  if (run == 201168 && (lumi >= 483 && lumi <= 558))
    return true;
  if (run == 201168 && (lumi >= 560 && lumi <= 730))
    return true;
  if (run == 201164 && (lumi >= 1 && lumi <= 8))
    return true;
  if (run == 201164 && (lumi >= 10 && lumi <= 94))
    return true;
  if (run == 201164 && (lumi >= 96 && lumi <= 125))
    return true;
  if (run == 201164 && (lumi >= 128 && lumi <= 178))
    return true;
  if (run == 201164 && (lumi >= 180 && lumi <= 198))
    return true;
  if (run == 201164 && (lumi >= 200 && lumi <= 271))
    return true;
  if (run == 201164 && (lumi >= 274 && lumi <= 416))
    return true;
  if (run == 201164 && (lumi >= 418 && lumi <= 418))
    return true;
  if (run == 195775 && (lumi >= 1 && lumi <= 57))
    return true;
  if (run == 195775 && (lumi >= 60 && lumi <= 100))
    return true;
  if (run == 195775 && (lumi >= 103 && lumi <= 170))
    return true;
  if (run == 195774 && (lumi >= 1 && lumi <= 13))
    return true;
  if (run == 195774 && (lumi >= 16 && lumi <= 137))
    return true;
  if (run == 195774 && (lumi >= 139 && lumi <= 151))
    return true;
  if (run == 195774 && (lumi >= 154 && lumi <= 162))
    return true;
  if (run == 195774 && (lumi >= 164 && lumi <= 256))
    return true;
  if (run == 195774 && (lumi >= 258 && lumi <= 276))
    return true;
  if (run == 195774 && (lumi >= 279 && lumi <= 362))
    return true;
  if (run == 195774 && (lumi >= 365 && lumi <= 466))
    return true;
  if (run == 195774 && (lumi >= 469 && lumi <= 618))
    return true;
  if (run == 195774 && (lumi >= 620 && lumi <= 649))
    return true;
  if (run == 195774 && (lumi >= 651 && lumi <= 830))
    return true;
  if (run == 195776 && (lumi >= 1 && lumi <= 63))
    return true;
  if (run == 195776 && (lumi >= 66 && lumi <= 283))
    return true;
  if (run == 195776 && (lumi >= 286 && lumi <= 337))
    return true;
  if (run == 195776 && (lumi >= 340 && lumi <= 399))
    return true;
  if (run == 195776 && (lumi >= 401 && lumi <= 409))
    return true;
  if (run == 195776 && (lumi >= 411 && lumi <= 477))
    return true;
  if (run == 194825 && (lumi >= 72 && lumi <= 117))
    return true;
  if (run == 194825 && (lumi >= 120 && lumi <= 221))
    return true;
  if (run == 206066 && (lumi >= 89 && lumi <= 146))
    return true;
  if (run == 206745 && (lumi >= 1 && lumi <= 81))
    return true;
  if (run == 206745 && (lumi >= 84 && lumi <= 199))
    return true;
  if (run == 206745 && (lumi >= 202 && lumi <= 224))
    return true;
  if (run == 206745 && (lumi >= 227 && lumi <= 237))
    return true;
  if (run == 206745 && (lumi >= 240 && lumi <= 304))
    return true;
  if (run == 206745 && (lumi >= 306 && lumi <= 318))
    return true;
  if (run == 206745 && (lumi >= 321 && lumi <= 720))
    return true;
  if (run == 206745 && (lumi >= 723 && lumi <= 796))
    return true;
  if (run == 206745 && (lumi >= 799 && lumi <= 894))
    return true;
  if (run == 206745 && (lumi >= 897 && lumi <= 944))
    return true;
  if (run == 206745 && (lumi >= 946 && lumi <= 1106))
    return true;
  if (run == 206745 && (lumi >= 1108 && lumi <= 1524))
    return true;
  if (run == 206745 && (lumi >= 1527 && lumi <= 1862))
    return true;
  if (run == 206745 && (lumi >= 1988 && lumi <= 1996))
    return true;
  if (run == 195147 && (lumi >= 132 && lumi <= 282))
    return true;
  if (run == 195147 && (lumi >= 285 && lumi <= 294))
    return true;
  if (run == 195147 && (lumi >= 297 && lumi <= 331))
    return true;
  if (run == 195147 && (lumi >= 334 && lumi <= 363))
    return true;
  if (run == 195147 && (lumi >= 366 && lumi <= 442))
    return true;
  if (run == 195147 && (lumi >= 445 && lumi <= 536))
    return true;
  if (run == 195147 && (lumi >= 539 && lumi <= 559))
    return true;
  if (run == 200369 && (lumi >= 1 && lumi <= 5))
    return true;
  if (run == 200369 && (lumi >= 8 && lumi <= 61))
    return true;
  if (run == 200369 && (lumi >= 64 && lumi <= 360))
    return true;
  if (run == 200369 && (lumi >= 363 && lumi <= 439))
    return true;
  if (run == 200369 && (lumi >= 441 && lumi <= 578))
    return true;
  if (run == 200369 && (lumi >= 580 && lumi <= 603))
    return true;
  if (run == 200369 && (lumi >= 606 && lumi <= 684))
    return true;
  if (run == 200369 && (lumi >= 686 && lumi <= 686))
    return true;
  if (run == 200368 && (lumi >= 72 && lumi <= 180))
    return true;
  if (run == 207488 && (lumi >= 1 && lumi <= 63))
    return true;
  if (run == 207488 && (lumi >= 66 && lumi <= 92))
    return true;
  if (run == 207488 && (lumi >= 95 && lumi <= 113))
    return true;
  if (run == 207488 && (lumi >= 116 && lumi <= 198))
    return true;
  if (run == 207488 && (lumi >= 200 && lumi <= 250))
    return true;
  if (run == 207488 && (lumi >= 252 && lumi <= 288))
    return true;
  if (run == 207488 && (lumi >= 291 && lumi <= 365))
    return true;
  if (run == 207488 && (lumi >= 368 && lumi <= 377))
    return true;
  if (run == 207488 && (lumi >= 379 && lumi <= 440))
    return true;
  if (run == 194439 && (lumi >= 46 && lumi <= 77))
    return true;
  if (run == 194439 && (lumi >= 79 && lumi <= 106))
    return true;
  if (run == 202389 && (lumi >= 81 && lumi <= 182))
    return true;
  if (run == 202389 && (lumi >= 185 && lumi <= 190))
    return true;
  if (run == 202389 && (lumi >= 192 && lumi <= 199))
    return true;
  if (run == 191859 && (lumi >= 1 && lumi <= 28))
    return true;
  if (run == 191859 && (lumi >= 31 && lumi <= 126))
    return true;
  if (run == 207518 && (lumi >= 1 && lumi <= 59))
    return true;
  if (run == 207518 && (lumi >= 61 && lumi <= 83))
    return true;
  if (run == 201554 && (lumi >= 70 && lumi <= 86))
    return true;
  if (run == 201554 && (lumi >= 88 && lumi <= 114))
    return true;
  if (run == 201554 && (lumi >= 116 && lumi <= 126))
    return true;
  if (run == 207517 && (lumi >= 1 && lumi <= 12))
    return true;
  if (run == 207517 && (lumi >= 15 && lumi <= 57))
    return true;
  if (run == 199011 && (lumi >= 1 && lumi <= 11))
    return true;
  if (run == 199011 && (lumi >= 13 && lumi <= 24))
    return true;
  if (run == 191856 && (lumi >= 1 && lumi <= 133))
    return true;
  if (run == 204555 && (lumi >= 1 && lumi <= 329))
    return true;
  if (run == 204555 && (lumi >= 331 && lumi <= 334))
    return true;
  if (run == 205086 && (lumi >= 95 && lumi <= 149))
    return true;
  if (run == 196197 && (lumi >= 58 && lumi <= 122))
    return true;
  if (run == 196197 && (lumi >= 125 && lumi <= 179))
    return true;
  if (run == 196197 && (lumi >= 181 && lumi <= 311))
    return true;
  if (run == 196197 && (lumi >= 313 && lumi <= 516))
    return true;
  if (run == 196197 && (lumi >= 519 && lumi <= 562))
    return true;
  if (run == 200992 && (lumi >= 1 && lumi <= 405))
    return true;
  if (run == 200992 && (lumi >= 408 && lumi <= 434))
    return true;
  if (run == 200992 && (lumi >= 436 && lumi <= 581))
    return true;
  if (run == 196199 && (lumi >= 1 && lumi <= 33))
    return true;
  if (run == 196199 && (lumi >= 36 && lumi <= 83))
    return true;
  if (run == 196199 && (lumi >= 86 && lumi <= 118))
    return true;
  if (run == 196199 && (lumi >= 121 && lumi <= 147))
    return true;
  if (run == 196199 && (lumi >= 150 && lumi <= 237))
    return true;
  if (run == 196199 && (lumi >= 239 && lumi <= 285))
    return true;
  if (run == 196199 && (lumi >= 287 && lumi <= 534))
    return true;
  if (run == 200990 && (lumi >= 75 && lumi <= 143))
    return true;
  if (run == 193334 && (lumi >= 29 && lumi <= 172))
    return true;
  if (run == 198269 && (lumi >= 3 && lumi <= 198))
    return true;
  if (run == 194075 && (lumi >= 48 && lumi <= 101))
    return true;
  if (run == 194075 && (lumi >= 103 && lumi <= 103))
    return true;
  if (run == 194075 && (lumi >= 105 && lumi <= 107))
    return true;
  if (run == 194075 && (lumi >= 109 && lumi <= 109))
    return true;
  if (run == 194075 && (lumi >= 111 && lumi <= 111))
    return true;
  if (run == 200519 && (lumi >= 1 && lumi <= 111))
    return true;
  if (run == 200519 && (lumi >= 114 && lumi <= 126))
    return true;
  if (run == 200519 && (lumi >= 129 && lumi <= 136))
    return true;
  if (run == 200519 && (lumi >= 138 && lumi <= 224))
    return true;
  if (run == 200519 && (lumi >= 227 && lumi <= 258))
    return true;
  if (run == 200519 && (lumi >= 261 && lumi <= 350))
    return true;
  if (run == 200519 && (lumi >= 353 && lumi <= 611))
    return true;
  if (run == 200519 && (lumi >= 613 && lumi <= 747))
    return true;
  if (run == 203742 && (lumi >= 1 && lumi <= 29))
    return true;
  if (run == 200515 && (lumi >= 97 && lumi <= 183))
    return true;
  if (run == 194643 && (lumi >= 1 && lumi <= 287))
    return true;
  if (run == 208487 && (lumi >= 2 && lumi <= 26))
    return true;
  if (run == 208487 && (lumi >= 29 && lumi <= 159))
    return true;
  if (run == 208487 && (lumi >= 161 && lumi <= 307))
    return true;
  if (run == 208487 && (lumi >= 309 && lumi <= 459))
    return true;
  if (run == 208487 && (lumi >= 462 && lumi <= 476))
    return true;
  if (run == 208487 && (lumi >= 479 && lumi <= 621))
    return true;
  if (run == 199751 && (lumi >= 103 && lumi <= 119))
    return true;
  if (run == 199751 && (lumi >= 121 && lumi <= 127))
    return true;
  if (run == 194644 && (lumi >= 1 && lumi <= 168))
    return true;
  if (run == 194644 && (lumi >= 171 && lumi <= 181))
    return true;
  if (run == 194644 && (lumi >= 184 && lumi <= 185))
    return true;
  if (run == 194644 && (lumi >= 187 && lumi <= 319))
    return true;
  if (run == 194644 && (lumi >= 321 && lumi <= 421))
    return true;
  if (run == 191411 && (lumi >= 1 && lumi <= 23))
    return true;
  if (run == 206246 && (lumi >= 1 && lumi <= 14))
    return true;
  if (run == 206246 && (lumi >= 16 && lumi <= 237))
    return true;
  if (run == 206246 && (lumi >= 240 && lumi <= 285))
    return true;
  if (run == 206246 && (lumi >= 288 && lumi <= 407))
    return true;
  if (run == 206246 && (lumi >= 412 && lumi <= 676))
    return true;
  if (run == 206246 && (lumi >= 678 && lumi <= 704))
    return true;
  if (run == 206246 && (lumi >= 706 && lumi <= 785))
    return true;
  if (run == 206246 && (lumi >= 787 && lumi <= 962))
    return true;
  if (run == 206246 && (lumi >= 965 && lumi <= 997))
    return true;
  if (run == 206246 && (lumi >= 1000 && lumi <= 1198))
    return true;
  if (run == 206246 && (lumi >= 1201 && lumi <= 1290))
    return true;
  if (run == 204599 && (lumi >= 73 && lumi <= 83))
    return true;
  if (run == 204599 && (lumi >= 85 && lumi <= 94))
    return true;
  if (run == 204599 && (lumi >= 97 && lumi <= 121))
    return true;
  if (run == 204599 && (lumi >= 124 && lumi <= 125))
    return true;
  if (run == 204599 && (lumi >= 128 && lumi <= 173))
    return true;
  if (run == 204599 && (lumi >= 175 && lumi <= 240))
    return true;
  if (run == 204599 && (lumi >= 243 && lumi <= 245))
    return true;
  if (run == 204599 && (lumi >= 248 && lumi <= 264))
    return true;
  if (run == 204599 && (lumi >= 266 && lumi <= 292))
    return true;
  if (run == 204599 && (lumi >= 294 && lumi <= 334))
    return true;
  if (run == 206243 && (lumi >= 62 && lumi <= 169))
    return true;
  if (run == 206243 && (lumi >= 172 && lumi <= 196))
    return true;
  if (run == 206243 && (lumi >= 199 && lumi <= 354))
    return true;
  if (run == 206243 && (lumi >= 357 && lumi <= 433))
    return true;
  if (run == 206243 && (lumi >= 435 && lumi <= 448))
    return true;
  if (run == 206243 && (lumi >= 451 && lumi <= 533))
    return true;
  if (run == 206243 && (lumi >= 536 && lumi <= 554))
    return true;
  if (run == 206243 && (lumi >= 557 && lumi <= 723))
    return true;
  if (run == 206243 && (lumi >= 726 && lumi <= 905))
    return true;
  if (run == 205344 && (lumi >= 1 && lumi <= 83))
    return true;
  if (run == 205344 && (lumi >= 86 && lumi <= 104))
    return true;
  if (run == 205344 && (lumi >= 106 && lumi <= 359))
    return true;
  if (run == 205344 && (lumi >= 362 && lumi <= 431))
    return true;
  if (run == 205344 && (lumi >= 433 && lumi <= 949))
    return true;
  if (run == 205344 && (lumi >= 951 && lumi <= 967))
    return true;
  if (run == 205344 && (lumi >= 969 && lumi <= 1127))
    return true;
  if (run == 205344 && (lumi >= 1129 && lumi <= 1346))
    return true;
  if (run == 205344 && (lumi >= 1348 && lumi <= 1586))
    return true;
  if (run == 200190 && (lumi >= 1 && lumi <= 4))
    return true;
  if (run == 200190 && (lumi >= 6 && lumi <= 76))
    return true;
  if (run == 200190 && (lumi >= 79 && lumi <= 143))
    return true;
  if (run == 200190 && (lumi >= 146 && lumi <= 159))
    return true;
  if (run == 200190 && (lumi >= 162 && lumi <= 256))
    return true;
  if (run == 200190 && (lumi >= 258 && lumi <= 321))
    return true;
  if (run == 200190 && (lumi >= 324 && lumi <= 401))
    return true;
  if (run == 200190 && (lumi >= 403 && lumi <= 453))
    return true;
  if (run == 200190 && (lumi >= 456 && lumi <= 457))
    return true;
  if (run == 200190 && (lumi >= 460 && lumi <= 565))
    return true;
  if (run == 200190 && (lumi >= 567 && lumi <= 588))
    return true;
  if (run == 200190 && (lumi >= 591 && lumi <= 591))
    return true;
  if (run == 200190 && (lumi >= 593 && lumi <= 595))
    return true;
  if (run == 200190 && (lumi >= 597 && lumi <= 646))
    return true;
  if (run == 200190 && (lumi >= 649 && lumi <= 878))
    return true;
  if (run == 204601 && (lumi >= 1 && lumi <= 25))
    return true;
  if (run == 204601 && (lumi >= 28 && lumi <= 62))
    return true;
  if (run == 204601 && (lumi >= 65 && lumi <= 80))
    return true;
  if (run == 204601 && (lumi >= 83 && lumi <= 89))
    return true;
  if (run == 204601 && (lumi >= 92 && lumi <= 290))
    return true;
  if (run == 204601 && (lumi >= 292 && lumi <= 563))
    return true;
  if (run == 204601 && (lumi >= 565 && lumi <= 591))
    return true;
  if (run == 204601 && (lumi >= 593 && lumi <= 652))
    return true;
  if (run == 204601 && (lumi >= 655 && lumi <= 780))
    return true;
  if (run == 204601 && (lumi >= 783 && lumi <= 812))
    return true;
  if (run == 204601 && (lumi >= 814 && lumi <= 892))
    return true;
  if (run == 204601 && (lumi >= 894 && lumi <= 984))
    return true;
  if (run == 204601 && (lumi >= 986 && lumi <= 1003))
    return true;
  if (run == 204601 && (lumi >= 1006 && lumi <= 1038))
    return true;
  if (run == 204601 && (lumi >= 1040 && lumi <= 1088))
    return true;
  if (run == 204601 && (lumi >= 1091 && lumi <= 1102))
    return true;
  if (run == 204601 && (lumi >= 1105 && lumi <= 1161))
    return true;
  if (run == 204601 && (lumi >= 1164 && lumi <= 1250))
    return true;
  if (run == 194480 && (lumi >= 1 && lumi <= 32))
    return true;
  if (run == 194480 && (lumi >= 34 && lumi <= 205))
    return true;
  if (run == 194480 && (lumi >= 207 && lumi <= 375))
    return true;
  if (run == 194480 && (lumi >= 377 && lumi <= 387))
    return true;
  if (run == 194480 && (lumi >= 389 && lumi <= 759))
    return true;
  if (run == 194480 && (lumi >= 762 && lumi <= 956))
    return true;
  if (run == 194480 && (lumi >= 959 && lumi <= 1402))
    return true;
  if (run == 205833 && (lumi >= 84 && lumi <= 86))
    return true;
  if (run == 205833 && (lumi >= 89 && lumi <= 121))
    return true;
  if (run == 205833 && (lumi >= 123 && lumi <= 155))
    return true;
  if (run == 205833 && (lumi >= 157 && lumi <= 165))
    return true;
  if (run == 205833 && (lumi >= 167 && lumi <= 173))
    return true;
  if (run == 205833 && (lumi >= 176 && lumi <= 219))
    return true;
  if (run == 205833 && (lumi >= 221 && lumi <= 267))
    return true;
  if (run == 205833 && (lumi >= 270 && lumi <= 312))
    return true;
  if (run == 205833 && (lumi >= 315 && lumi <= 346))
    return true;
  if (run == 205833 && (lumi >= 350 && lumi <= 355))
    return true;
  if (run == 205833 && (lumi >= 360 && lumi <= 366))
    return true;
  if (run == 205834 && (lumi >= 1 && lumi <= 12))
    return true;
  if (run == 205834 && (lumi >= 14 && lumi <= 195))
    return true;
  if (run == 193621 && (lumi >= 60 && lumi <= 570))
    return true;
  if (run == 193621 && (lumi >= 573 && lumi <= 769))
    return true;
  if (run == 193621 && (lumi >= 772 && lumi <= 976))
    return true;
  if (run == 193621 && (lumi >= 979 && lumi <= 1053))
    return true;
  if (run == 193621 && (lumi >= 1056 && lumi <= 1137))
    return true;
  if (run == 193621 && (lumi >= 1139 && lumi <= 1193))
    return true;
  if (run == 193621 && (lumi >= 1195 && lumi <= 1371))
    return true;
  if (run == 193621 && (lumi >= 1373 && lumi <= 1654))
    return true;
  if (run == 191090 && (lumi >= 1 && lumi <= 55))
    return true;
  if (run == 205515 && (lumi >= 82 && lumi <= 201))
    return true;
  if (run == 205515 && (lumi >= 203 && lumi <= 216))
    return true;
  if (run == 205519 && (lumi >= 1 && lumi <= 47))
    return true;
  if (run == 205519 && (lumi >= 50 && lumi <= 172))
    return true;
  if (run == 205519 && (lumi >= 175 && lumi <= 367))
    return true;
  if (run == 205519 && (lumi >= 370 && lumi <= 386))
    return true;
  if (run == 205519 && (lumi >= 389 && lumi <= 472))
    return true;
  if (run == 194897 && (lumi >= 1 && lumi <= 6))
    return true;
  if (run == 194897 && (lumi >= 8 && lumi <= 78))
    return true;
  if (run == 194897 && (lumi >= 80 && lumi <= 96))
    return true;
  if (run == 194897 && (lumi >= 98 && lumi <= 102))
    return true;
  if (run == 194896 && (lumi >= 34 && lumi <= 55))
    return true;
  if (run == 194896 && (lumi >= 58 && lumi <= 79))
    return true;
  if (run == 194896 && (lumi >= 82 && lumi <= 103))
    return true;
  if (run == 201602 && (lumi >= 76 && lumi <= 81))
    return true;
  if (run == 201602 && (lumi >= 83 && lumi <= 194))
    return true;
  if (run == 201602 && (lumi >= 196 && lumi <= 494))
    return true;
  if (run == 201602 && (lumi >= 496 && lumi <= 614))
    return true;
  if (run == 201602 && (lumi >= 617 && lumi <= 635))
    return true;
  if (run == 207515 && (lumi >= 79 && lumi <= 109))
    return true;
  if (run == 207515 && (lumi >= 112 && lumi <= 132))
    return true;
  if (run == 207515 && (lumi >= 134 && lumi <= 208))
    return true;
  if (run == 207515 && (lumi >= 211 && lumi <= 225))
    return true;
  if (run == 207515 && (lumi >= 228 && lumi <= 320))
    return true;
  if (run == 207515 && (lumi >= 322 && lumi <= 381))
    return true;
  if (run == 207515 && (lumi >= 383 && lumi <= 498))
    return true;
  if (run == 207515 && (lumi >= 500 && lumi <= 730))
    return true;
  if (run == 207515 && (lumi >= 733 && lumi <= 849))
    return true;
  if (run == 207515 && (lumi >= 851 && lumi <= 954))
    return true;
  if (run == 207515 && (lumi >= 957 && lumi <= 994))
    return true;
  if (run == 207515 && (lumi >= 997 && lumi <= 1052))
    return true;
  if (run == 207515 && (lumi >= 1055 && lumi <= 1143))
    return true;
  if (run == 207515 && (lumi >= 1145 && lumi <= 1211))
    return true;
  if (run == 195923 && (lumi >= 1 && lumi <= 14))
    return true;
  if (run == 195925 && (lumi >= 1 && lumi <= 12))
    return true;
  if (run == 195926 && (lumi >= 1 && lumi <= 1))
    return true;
  if (run == 195926 && (lumi >= 3 && lumi <= 19))
    return true;
  if (run == 195926 && (lumi >= 21 && lumi <= 34))
    return true;
  if (run == 208351 && (lumi >= 83 && lumi <= 97))
    return true;
  if (run == 208351 && (lumi >= 100 && lumi <= 356))
    return true;
  if (run == 208351 && (lumi >= 359 && lumi <= 367))
    return true;
  if (run == 208351 && (lumi >= 369 && lumi <= 369))
    return true;
  if (run == 208352 && (lumi >= 1 && lumi <= 15))
    return true;
  if (run == 208352 && (lumi >= 17 && lumi <= 17))
    return true;
  if (run == 208352 && (lumi >= 19 && lumi <= 19))
    return true;
  if (run == 208353 && (lumi >= 1 && lumi <= 76))
    return true;
  if (run == 208353 && (lumi >= 78 && lumi <= 269))
    return true;
  if (run == 208353 && (lumi >= 271 && lumi <= 348))
    return true;
  if (run == 208357 && (lumi >= 1 && lumi <= 70))
    return true;
  if (run == 208357 && (lumi >= 73 && lumi <= 507))
    return true;
  if (run == 199832 && (lumi >= 58 && lumi <= 62))
    return true;
  if (run == 199832 && (lumi >= 65 && lumi <= 118))
    return true;
  if (run == 199832 && (lumi >= 121 && lumi <= 139))
    return true;
  if (run == 199832 && (lumi >= 142 && lumi <= 286))
    return true;
  if (run == 199833 && (lumi >= 1 && lumi <= 13))
    return true;
  if (run == 199833 && (lumi >= 16 && lumi <= 103))
    return true;
  if (run == 199833 && (lumi >= 105 && lumi <= 250))
    return true;
  if (run == 199833 && (lumi >= 253 && lumi <= 493))
    return true;
  if (run == 199833 && (lumi >= 496 && lumi <= 794))
    return true;
  if (run == 199833 && (lumi >= 797 && lumi <= 1032))
    return true;
  if (run == 199833 && (lumi >= 1034 && lumi <= 1185))
    return true;
  if (run == 199833 && (lumi >= 1188 && lumi <= 1239))
    return true;
  if (run == 207469 && (lumi >= 1 && lumi <= 31))
    return true;
  if (run == 207469 && (lumi >= 34 && lumi <= 45))
    return true;
  if (run == 199834 && (lumi >= 1 && lumi <= 9))
    return true;
  if (run == 199834 && (lumi >= 11 && lumi <= 11))
    return true;
  if (run == 199834 && (lumi >= 14 && lumi <= 18))
    return true;
  if (run == 199834 && (lumi >= 21 && lumi <= 54))
    return true;
  if (run == 199834 && (lumi >= 56 && lumi <= 57))
    return true;
  if (run == 199834 && (lumi >= 62 && lumi <= 65))
    return true;
  if (run == 199834 && (lumi >= 69 && lumi <= 284))
    return true;
  if (run == 199834 && (lumi >= 286 && lumi <= 503))
    return true;
  if (run == 199834 && (lumi >= 505 && lumi <= 942))
    return true;
  if (run == 194790 && (lumi >= 1 && lumi <= 45))
    return true;
  if (run == 204114 && (lumi >= 1 && lumi <= 358))
    return true;
  if (run == 206391 && (lumi >= 1 && lumi <= 55))
    return true;
  if (run == 206391 && (lumi >= 57 && lumi <= 91))
    return true;
  if (run == 207269 && (lumi >= 80 && lumi <= 394))
    return true;
  if (run == 207269 && (lumi >= 397 && lumi <= 436))
    return true;
  if (run == 207269 && (lumi >= 439 && lumi <= 463))
    return true;
  if (run == 207269 && (lumi >= 466 && lumi <= 551))
    return true;
  if (run == 207269 && (lumi >= 568 && lumi <= 577))
    return true;
  if (run == 204113 && (lumi >= 82 && lumi <= 96))
    return true;
  if (run == 204113 && (lumi >= 98 && lumi <= 102))
    return true;
  if (run == 204113 && (lumi >= 105 && lumi <= 127))
    return true;
  if (run == 204113 && (lumi >= 129 && lumi <= 191))
    return true;
  if (run == 204113 && (lumi >= 194 && lumi <= 258))
    return true;
  if (run == 204113 && (lumi >= 261 && lumi <= 327))
    return true;
  if (run == 204113 && (lumi >= 329 && lumi <= 388))
    return true;
  if (run == 204113 && (lumi >= 390 && lumi <= 400))
    return true;
  if (run == 204113 && (lumi >= 402 && lumi <= 583))
    return true;
  if (run == 204113 && (lumi >= 585 && lumi <= 690))
    return true;
  if (run == 204238 && (lumi >= 23 && lumi <= 52))
    return true;
  if (run == 204238 && (lumi >= 55 && lumi <= 55))
    return true;
  if (run == 195633 && (lumi >= 40 && lumi <= 42))
    return true;
  if (run == 206897 && (lumi >= 1 && lumi <= 34))
    return true;
  if (run == 206897 && (lumi >= 38 && lumi <= 61))
    return true;
  if (run == 206897 && (lumi >= 63 && lumi <= 102))
    return true;
  if (run == 206897 && (lumi >= 109 && lumi <= 109))
    return true;
  if (run == 206897 && (lumi >= 111 && lumi <= 112))
    return true;
  if (run == 206897 && (lumi >= 114 && lumi <= 131))
    return true;
  if (run == 206897 && (lumi >= 133 && lumi <= 137))
    return true;
  if (run == 196249 && (lumi >= 63 && lumi <= 77))
    return true;
  if (run == 196249 && (lumi >= 80 && lumi <= 99))
    return true;
  if (run == 200991 && (lumi >= 1 && lumi <= 42))
    return true;
  if (run == 200991 && (lumi >= 44 && lumi <= 44))
    return true;
  if (run == 200991 && (lumi >= 47 && lumi <= 80))
    return true;
  if (run == 200991 && (lumi >= 83 && lumi <= 175))
    return true;
  if (run == 200991 && (lumi >= 178 && lumi <= 181))
    return true;
  if (run == 200991 && (lumi >= 184 && lumi <= 252))
    return true;
  if (run == 200991 && (lumi >= 255 && lumi <= 632))
    return true;
  if (run == 200991 && (lumi >= 635 && lumi <= 916))
    return true;
  if (run == 200991 && (lumi >= 918 && lumi <= 1017))
    return true;
  if (run == 200991 && (lumi >= 1019 && lumi <= 1048))
    return true;
  if (run == 194210 && (lumi >= 3 && lumi <= 195))
    return true;
  if (run == 194210 && (lumi >= 198 && lumi <= 217))
    return true;
  if (run == 194210 && (lumi >= 220 && lumi <= 359))
    return true;
  if (run == 194210 && (lumi >= 361 && lumi <= 555))
    return true;
  if (run == 208509 && (lumi >= 71 && lumi <= 232))
    return true;
  if (run == 200049 && (lumi >= 1 && lumi <= 177))
    return true;
  if (run == 198116 && (lumi >= 36 && lumi <= 52))
    return true;
  if (run == 198116 && (lumi >= 54 && lumi <= 55))
    return true;
  if (run == 198116 && (lumi >= 58 && lumi <= 96))
    return true;
  if (run == 198116 && (lumi >= 98 && lumi <= 112))
    return true;
  if (run == 203853 && (lumi >= 122 && lumi <= 222))
    return true;
  if (run == 190688 && (lumi >= 69 && lumi <= 249))
    return true;
  if (run == 206245 && (lumi >= 1 && lumi <= 62))
    return true;
  if (run == 207492 && (lumi >= 1 && lumi <= 20))
    return true;
  if (run == 207492 && (lumi >= 23 && lumi <= 298))
    return true;
  if (run == 206210 && (lumi >= 11 && lumi <= 25))
    return true;
  if (run == 206210 && (lumi >= 28 && lumi <= 275))
    return true;
  if (run == 206210 && (lumi >= 277 && lumi <= 298))
    return true;
  if (run == 206210 && (lumi >= 300 && lumi <= 383))
    return true;
  if (run == 206210 && (lumi >= 386 && lumi <= 466))
    return true;
  if (run == 198230 && (lumi >= 1 && lumi <= 33))
    return true;
  if (run == 198230 && (lumi >= 36 && lumi <= 57))
    return true;
  if (run == 198230 && (lumi >= 60 && lumi <= 235))
    return true;
  if (run == 198230 && (lumi >= 237 && lumi <= 324))
    return true;
  if (run == 198230 && (lumi >= 326 && lumi <= 388))
    return true;
  if (run == 198230 && (lumi >= 390 && lumi <= 459))
    return true;
  if (run == 198230 && (lumi >= 462 && lumi <= 625))
    return true;
  if (run == 198230 && (lumi >= 627 && lumi <= 651))
    return true;
  if (run == 198230 && (lumi >= 653 && lumi <= 805))
    return true;
  if (run == 198230 && (lumi >= 808 && lumi <= 811))
    return true;
  if (run == 198230 && (lumi >= 814 && lumi <= 948))
    return true;
  if (run == 198230 && (lumi >= 950 && lumi <= 1090))
    return true;
  if (run == 198230 && (lumi >= 1093 && lumi <= 1103))
    return true;
  if (run == 198230 && (lumi >= 1106 && lumi <= 1332))
    return true;
  if (run == 198230 && (lumi >= 1335 && lumi <= 1380))
    return true;
  if (run == 208341 && (lumi >= 1 && lumi <= 84))
    return true;
  if (run == 208341 && (lumi >= 87 && lumi <= 117))
    return true;
  if (run == 208341 && (lumi >= 120 && lumi <= 513))
    return true;
  if (run == 208341 && (lumi >= 515 && lumi <= 685))
    return true;
  if (run == 208341 && (lumi >= 688 && lumi <= 693))
    return true;
  if (run == 208341 && (lumi >= 695 && lumi <= 775))
    return true;
  if (run == 208341 && (lumi >= 777 && lumi <= 824))
    return true;
  if (run == 191810 && (lumi >= 15 && lumi <= 15))
    return true;
  if (run == 191810 && (lumi >= 22 && lumi <= 49))
    return true;
  if (run == 191810 && (lumi >= 52 && lumi <= 92))
    return true;
  if (run == 208428 && (lumi >= 1 && lumi <= 58))
    return true;
  if (run == 208428 && (lumi >= 61 && lumi <= 68))
    return true;
  if (run == 208428 && (lumi >= 70 && lumi <= 156))
    return true;
  if (run == 208428 && (lumi >= 159 && lumi <= 227))
    return true;
  if (run == 196019 && (lumi >= 54 && lumi <= 68))
    return true;
  if (run == 198941 && (lumi >= 102 && lumi <= 189))
    return true;
  if (run == 198941 && (lumi >= 191 && lumi <= 220))
    return true;
  if (run == 198941 && (lumi >= 222 && lumi <= 241))
    return true;
  if (run == 198941 && (lumi >= 243 && lumi <= 249))
    return true;
  if (run == 198941 && (lumi >= 252 && lumi <= 284))
    return true;
  if (run == 194076 && (lumi >= 1 && lumi <= 9))
    return true;
  if (run == 194076 && (lumi >= 11 && lumi <= 55))
    return true;
  if (run == 194076 && (lumi >= 58 && lumi <= 163))
    return true;
  if (run == 194076 && (lumi >= 165 && lumi <= 228))
    return true;
  if (run == 194076 && (lumi >= 230 && lumi <= 264))
    return true;
  if (run == 194076 && (lumi >= 267 && lumi <= 507))
    return true;
  if (run == 194076 && (lumi >= 509 && lumi <= 527))
    return true;
  if (run == 194076 && (lumi >= 530 && lumi <= 538))
    return true;
  if (run == 194076 && (lumi >= 541 && lumi <= 562))
    return true;
  if (run == 194076 && (lumi >= 565 && lumi <= 748))
    return true;
  if (run == 195115 && (lumi >= 1 && lumi <= 7))
    return true;
  if (run == 195115 && (lumi >= 10 && lumi <= 22))
    return true;
  if (run == 195114 && (lumi >= 1 && lumi <= 69))
    return true;
  if (run == 195114 && (lumi >= 72 && lumi <= 103))
    return true;
  if (run == 202299 && (lumi >= 68 && lumi <= 84))
    return true;
  if (run == 202299 && (lumi >= 87 && lumi <= 141))
    return true;
  if (run == 202299 && (lumi >= 143 && lumi <= 193))
    return true;
  if (run == 202299 && (lumi >= 196 && lumi <= 358))
    return true;
  if (run == 202299 && (lumi >= 361 && lumi <= 379))
    return true;
  if (run == 202299 && (lumi >= 382 && lumi <= 414))
    return true;
  if (run == 202299 && (lumi >= 416 && lumi <= 452))
    return true;
  if (run == 202299 && (lumi >= 455 && lumi <= 555))
    return true;
  if (run == 198346 && (lumi >= 44 && lumi <= 47))
    return true;
  if (run == 195113 && (lumi >= 1 && lumi <= 209))
    return true;
  if (run == 195113 && (lumi >= 212 && lumi <= 388))
    return true;
  if (run == 195113 && (lumi >= 391 && lumi <= 403))
    return true;
  if (run == 195113 && (lumi >= 406 && lumi <= 419))
    return true;
  if (run == 195113 && (lumi >= 422 && lumi <= 492))
    return true;
  if (run == 195113 && (lumi >= 495 && lumi <= 579))
    return true;
  if (run == 195112 && (lumi >= 1 && lumi <= 12))
    return true;
  if (run == 195112 && (lumi >= 15 && lumi <= 26))
    return true;
  if (run == 202075 && (lumi >= 1 && lumi <= 18))
    return true;
  if (run == 202075 && (lumi >= 21 && lumi <= 187))
    return true;
  if (run == 202075 && (lumi >= 189 && lumi <= 214))
    return true;
  if (run == 202075 && (lumi >= 217 && lumi <= 247))
    return true;
  if (run == 202075 && (lumi >= 250 && lumi <= 342))
    return true;
  if (run == 202075 && (lumi >= 345 && lumi <= 406))
    return true;
  if (run == 202075 && (lumi >= 409 && lumi <= 497))
    return true;
  if (run == 202075 && (lumi >= 500 && lumi <= 537))
    return true;
  if (run == 202075 && (lumi >= 539 && lumi <= 539))
    return true;
  if (run == 202075 && (lumi >= 542 && lumi <= 560))
    return true;
  if (run == 202075 && (lumi >= 562 && lumi <= 615))
    return true;
  if (run == 202075 && (lumi >= 618 && lumi <= 628))
    return true;
  if (run == 202074 && (lumi >= 66 && lumi <= 174))
    return true;
  if (run == 202178 && (lumi >= 67 && lumi <= 78))
    return true;
  if (run == 202178 && (lumi >= 80 && lumi <= 88))
    return true;
  if (run == 202178 && (lumi >= 91 && lumi <= 177))
    return true;
  if (run == 202178 && (lumi >= 180 && lumi <= 186))
    return true;
  if (run == 202178 && (lumi >= 188 && lumi <= 337))
    return true;
  if (run == 202178 && (lumi >= 340 && lumi <= 377))
    return true;
  if (run == 202178 && (lumi >= 379 && lumi <= 425))
    return true;
  if (run == 202178 && (lumi >= 428 && lumi <= 475))
    return true;
  if (run == 202178 && (lumi >= 478 && lumi <= 548))
    return true;
  if (run == 202178 && (lumi >= 551 && lumi <= 717))
    return true;
  if (run == 202178 && (lumi >= 720 && lumi <= 965))
    return true;
  if (run == 202178 && (lumi >= 967 && lumi <= 1444))
    return true;
  if (run == 202178 && (lumi >= 1447 && lumi <= 1505))
    return true;
  if (run == 202178 && (lumi >= 1508 && lumi <= 1519))
    return true;
  if (run == 202178 && (lumi >= 1522 && lumi <= 1555))
    return true;
  if (run == 202272 && (lumi >= 76 && lumi <= 112))
    return true;
  if (run == 202272 && (lumi >= 115 && lumi <= 141))
    return true;
  if (run == 202272 && (lumi >= 144 && lumi <= 185))
    return true;
  if (run == 202272 && (lumi >= 188 && lumi <= 205))
    return true;
  if (run == 202272 && (lumi >= 208 && lumi <= 305))
    return true;
  if (run == 202272 && (lumi >= 307 && lumi <= 313))
    return true;
  if (run == 202272 && (lumi >= 315 && lumi <= 371))
    return true;
  if (run == 202272 && (lumi >= 436 && lumi <= 480))
    return true;
  if (run == 202272 && (lumi >= 483 && lumi <= 555))
    return true;
  if (run == 202272 && (lumi >= 558 && lumi <= 577))
    return true;
  if (run == 202272 && (lumi >= 579 && lumi <= 683))
    return true;
  if (run == 202272 && (lumi >= 686 && lumi <= 705))
    return true;
  if (run == 202272 && (lumi >= 707 && lumi <= 740))
    return true;
  if (run == 202272 && (lumi >= 742 && lumi <= 890))
    return true;
  if (run == 202272 && (lumi >= 937 && lumi <= 1295))
    return true;
  if (run == 202272 && (lumi >= 1299 && lumi <= 1481))
    return true;
  if (run == 193541 && (lumi >= 77 && lumi <= 101))
    return true;
  if (run == 193541 && (lumi >= 103 && lumi <= 413))
    return true;
  if (run == 193541 && (lumi >= 416 && lumi <= 575))
    return true;
  if (run == 193541 && (lumi >= 578 && lumi <= 619))
    return true;
  if (run == 202328 && (lumi >= 46 && lumi <= 89))
    return true;
  if (run == 202328 && (lumi >= 92 && lumi <= 156))
    return true;
  if (run == 202328 && (lumi >= 158 && lumi <= 276))
    return true;
  if (run == 202328 && (lumi >= 278 && lumi <= 291))
    return true;
  if (run == 202328 && (lumi >= 294 && lumi <= 434))
    return true;
  if (run == 202328 && (lumi >= 437 && lumi <= 460))
    return true;
  if (run == 202328 && (lumi >= 463 && lumi <= 586))
    return true;
  if (run == 202328 && (lumi >= 588 && lumi <= 610))
    return true;
  if (run == 202328 && (lumi >= 612 && lumi <= 614))
    return true;
  if (run == 200174 && (lumi >= 81 && lumi <= 84))
    return true;
  if (run == 200177 && (lumi >= 1 && lumi <= 56))
    return true;
  if (run == 200178 && (lumi >= 1 && lumi <= 38))
    return true;
  if (run == 200075 && (lumi >= 76 && lumi <= 139))
    return true;
  if (run == 200075 && (lumi >= 142 && lumi <= 232))
    return true;
  if (run == 200075 && (lumi >= 256 && lumi <= 326))
    return true;
  if (run == 200075 && (lumi >= 329 && lumi <= 422))
    return true;
  if (run == 200075 && (lumi >= 425 && lumi <= 431))
    return true;
  if (run == 200075 && (lumi >= 434 && lumi <= 500))
    return true;
  if (run == 200075 && (lumi >= 502 && lumi <= 605))
    return true;
  if (run == 195963 && (lumi >= 54 && lumi <= 58))
    return true;
  if (run == 206199 && (lumi >= 1 && lumi <= 75))
    return true;
  if (run == 206199 && (lumi >= 77 && lumi <= 82))
    return true;
  if (run == 206199 && (lumi >= 85 && lumi <= 114))
    return true;
  if (run == 201816 && (lumi >= 1 && lumi <= 72))
    return true;
  if (run == 201816 && (lumi >= 74 && lumi <= 105))
    return true;
  if (run == 201816 && (lumi >= 107 && lumi <= 157))
    return true;
  if (run == 201817 && (lumi >= 1 && lumi <= 274))
    return true;
  if (run == 201818 && (lumi >= 1 && lumi <= 1))
    return true;
  if (run == 201819 && (lumi >= 1 && lumi <= 94))
    return true;
  if (run == 201819 && (lumi >= 96 && lumi <= 241))
    return true;
  if (run == 196349 && (lumi >= 65 && lumi <= 84))
    return true;
  if (run == 196349 && (lumi >= 86 && lumi <= 154))
    return true;
  if (run == 196349 && (lumi >= 157 && lumi <= 244))
    return true;
  if (run == 196349 && (lumi >= 246 && lumi <= 258))
    return true;
  if (run == 205908 && (lumi >= 68 && lumi <= 200))
    return true;
  if (run == 205908 && (lumi >= 202 && lumi <= 209))
    return true;
  if (run == 194619 && (lumi >= 31 && lumi <= 110))
    return true;
  if (run == 206605 && (lumi >= 1 && lumi <= 36))
    return true;
  if (run == 206605 && (lumi >= 39 && lumi <= 78))
    return true;
  if (run == 207221 && (lumi >= 1 && lumi <= 102))
    return true;
  if (run == 207220 && (lumi >= 1 && lumi <= 160))
    return true;
  if (run == 200491 && (lumi >= 87 && lumi <= 107))
    return true;
  if (run == 200491 && (lumi >= 110 && lumi <= 149))
    return true;
  if (run == 200491 && (lumi >= 152 && lumi <= 157))
    return true;
  if (run == 200491 && (lumi >= 160 && lumi <= 197))
    return true;
  if (run == 200491 && (lumi >= 199 && lumi <= 237))
    return true;
  if (run == 200491 && (lumi >= 240 && lumi <= 270))
    return true;
  if (run == 200491 && (lumi >= 273 && lumi <= 273))
    return true;
  if (run == 200491 && (lumi >= 276 && lumi <= 334))
    return true;
  if (run == 200491 && (lumi >= 336 && lumi <= 360))
    return true;
  if (run == 200491 && (lumi >= 363 && lumi <= 419))
    return true;
  if (run == 207222 && (lumi >= 1 && lumi <= 17))
    return true;
  if (run == 207222 && (lumi >= 20 && lumi <= 289))
    return true;
  if (run == 202209 && (lumi >= 1 && lumi <= 48))
    return true;
  if (run == 202209 && (lumi >= 51 && lumi <= 142))
    return true;
  if (run == 205311 && (lumi >= 1 && lumi <= 85))
    return true;
  if (run == 205311 && (lumi >= 88 && lumi <= 92))
    return true;
  if (run == 205311 && (lumi >= 95 && lumi <= 183))
    return true;
  if (run == 205311 && (lumi >= 186 && lumi <= 395))
    return true;
  if (run == 205311 && (lumi >= 397 && lumi <= 592))
    return true;
  if (run == 205311 && (lumi >= 595 && lumi <= 910))
    return true;
  if (run == 205311 && (lumi >= 913 && lumi <= 1260))
    return true;
  if (run == 205310 && (lumi >= 76 && lumi <= 306))
    return true;
  if (run == 205310 && (lumi >= 309 && lumi <= 313))
    return true;
  if (run == 205310 && (lumi >= 316 && lumi <= 316))
    return true;
  if (run == 205310 && (lumi >= 319 && lumi <= 321))
    return true;
  if (run == 205310 && (lumi >= 324 && lumi <= 457))
    return true;
  if (run == 205310 && (lumi >= 460 && lumi <= 559))
    return true;
  if (run == 206859 && (lumi >= 79 && lumi <= 210))
    return true;
  if (run == 206859 && (lumi >= 212 && lumi <= 258))
    return true;
  if (run == 206859 && (lumi >= 260 && lumi <= 323))
    return true;
  if (run == 206859 && (lumi >= 325 && lumi <= 356))
    return true;
  if (run == 206859 && (lumi >= 359 && lumi <= 609))
    return true;
  if (run == 206859 && (lumi >= 612 && lumi <= 681))
    return true;
  if (run == 206859 && (lumi >= 684 && lumi <= 732))
    return true;
  if (run == 206859 && (lumi >= 734 && lumi <= 768))
    return true;
  if (run == 206859 && (lumi >= 771 && lumi <= 808))
    return true;
  if (run == 206859 && (lumi >= 811 && lumi <= 827))
    return true;
  if (run == 206859 && (lumi >= 830 && lumi <= 848))
    return true;
  if (run == 194314 && (lumi >= 52 && lumi <= 130))
    return true;
  if (run == 194314 && (lumi >= 133 && lumi <= 300))
    return true;
  if (run == 194315 && (lumi >= 1 && lumi <= 10))
    return true;
  if (run == 194315 && (lumi >= 13 && lumi <= 314))
    return true;
  if (run == 194315 && (lumi >= 317 && lumi <= 428))
    return true;
  if (run == 194315 && (lumi >= 431 && lumi <= 452))
    return true;
  if (run == 194315 && (lumi >= 455 && lumi <= 467))
    return true;
  if (run == 194317 && (lumi >= 1 && lumi <= 20))
    return true;
  if (run == 200473 && (lumi >= 96 && lumi <= 157))
    return true;
  if (run == 200473 && (lumi >= 159 && lumi <= 224))
    return true;
  if (run == 200473 && (lumi >= 226 && lumi <= 304))
    return true;
  if (run == 200473 && (lumi >= 306 && lumi <= 469))
    return true;
  if (run == 200473 && (lumi >= 472 && lumi <= 524))
    return true;
  if (run == 200473 && (lumi >= 527 && lumi <= 542))
    return true;
  if (run == 200473 && (lumi >= 545 && lumi <= 619))
    return true;
  if (run == 200473 && (lumi >= 622 && lumi <= 688))
    return true;
  if (run == 200473 && (lumi >= 691 && lumi <= 730))
    return true;
  if (run == 200473 && (lumi >= 733 && lumi <= 738))
    return true;
  if (run == 200473 && (lumi >= 740 && lumi <= 1324))
    return true;
  if (run == 196200 && (lumi >= 1 && lumi <= 68))
    return true;
  if (run == 196202 && (lumi >= 3 && lumi <= 61))
    return true;
  if (run == 196202 && (lumi >= 64 && lumi <= 108))
    return true;
  if (run == 196203 && (lumi >= 1 && lumi <= 102))
    return true;
  if (run == 196203 && (lumi >= 107 && lumi <= 117))
    return true;
  if (run == 201678 && (lumi >= 1 && lumi <= 120))
    return true;
  if (run == 199864 && (lumi >= 1 && lumi <= 87))
    return true;
  if (run == 199864 && (lumi >= 89 && lumi <= 89))
    return true;
  if (run == 199864 && (lumi >= 92 && lumi <= 103))
    return true;
  if (run == 199864 && (lumi >= 106 && lumi <= 372))
    return true;
  if (run == 199864 && (lumi >= 374 && lumi <= 385))
    return true;
  if (run == 199864 && (lumi >= 388 && lumi <= 486))
    return true;
  if (run == 199867 && (lumi >= 1 && lumi <= 134))
    return true;
  if (run == 199867 && (lumi >= 136 && lumi <= 172))
    return true;
  if (run == 199867 && (lumi >= 174 && lumi <= 218))
    return true;
  if (run == 199867 && (lumi >= 221 && lumi <= 320))
    return true;
  if (run == 199862 && (lumi >= 59 && lumi <= 141))
    return true;
  if (run == 208541 && (lumi >= 1 && lumi <= 57))
    return true;
  if (run == 208541 && (lumi >= 59 && lumi <= 173))
    return true;
  if (run == 208541 && (lumi >= 175 && lumi <= 376))
    return true;
  if (run == 208541 && (lumi >= 378 && lumi <= 413))
    return true;
  if (run == 201671 && (lumi >= 1 && lumi <= 120))
    return true;
  if (run == 201671 && (lumi >= 122 && lumi <= 174))
    return true;
  if (run == 201671 && (lumi >= 177 && lumi <= 462))
    return true;
  if (run == 201671 && (lumi >= 464 && lumi <= 482))
    return true;
  if (run == 201671 && (lumi >= 485 && lumi <= 499))
    return true;
  if (run == 201671 && (lumi >= 501 && lumi <= 545))
    return true;
  if (run == 201671 && (lumi >= 547 && lumi <= 571))
    return true;
  if (run == 201671 && (lumi >= 574 && lumi <= 614))
    return true;
  if (run == 201671 && (lumi >= 617 && lumi <= 766))
    return true;
  if (run == 201671 && (lumi >= 768 && lumi <= 896))
    return true;
  if (run == 201671 && (lumi >= 899 && lumi <= 911))
    return true;
  if (run == 201671 && (lumi >= 914 && lumi <= 1007))
    return true;
  if (run == 203894 && (lumi >= 82 && lumi <= 272))
    return true;
  if (run == 203894 && (lumi >= 275 && lumi <= 477))
    return true;
  if (run == 203894 && (lumi >= 480 && lumi <= 902))
    return true;
  if (run == 203894 && (lumi >= 905 && lumi <= 1319))
    return true;
  if (run == 199868 && (lumi >= 1 && lumi <= 21))
    return true;
  if (run == 203777 && (lumi >= 103 && lumi <= 113))
    return true;
  if (run == 199319 && (lumi >= 1 && lumi <= 7))
    return true;
  if (run == 199319 && (lumi >= 9 && lumi <= 223))
    return true;
  if (run == 199319 && (lumi >= 226 && lumi <= 277))
    return true;
  if (run == 199319 && (lumi >= 280 && lumi <= 348))
    return true;
  if (run == 199319 && (lumi >= 351 && lumi <= 358))
    return true;
  if (run == 199319 && (lumi >= 360 && lumi <= 422))
    return true;
  if (run == 199319 && (lumi >= 424 && lumi <= 490))
    return true;
  if (run == 199319 && (lumi >= 492 && lumi <= 493))
    return true;
  if (run == 199319 && (lumi >= 496 && lumi <= 612))
    return true;
  if (run == 199319 && (lumi >= 615 && lumi <= 642))
    return true;
  if (run == 199319 && (lumi >= 645 && lumi <= 720))
    return true;
  if (run == 199319 && (lumi >= 723 && lumi <= 728))
    return true;
  if (run == 199319 && (lumi >= 730 && lumi <= 731))
    return true;
  if (run == 199319 && (lumi >= 734 && lumi <= 741))
    return true;
  if (run == 199319 && (lumi >= 744 && lumi <= 752))
    return true;
  if (run == 199319 && (lumi >= 754 && lumi <= 943))
    return true;
  if (run == 199319 && (lumi >= 945 && lumi <= 997))
    return true;
  if (run == 199703 && (lumi >= 1 && lumi <= 94))
    return true;
  if (run == 199703 && (lumi >= 97 && lumi <= 482))
    return true;
  if (run == 199703 && (lumi >= 485 && lumi <= 529))
    return true;
  if (run == 206512 && (lumi >= 91 && lumi <= 123))
    return true;
  if (run == 206512 && (lumi >= 125 && lumi <= 133))
    return true;
  if (run == 206512 && (lumi >= 136 && lumi <= 161))
    return true;
  if (run == 206512 && (lumi >= 163 && lumi <= 190))
    return true;
  if (run == 206512 && (lumi >= 193 && lumi <= 201))
    return true;
  if (run == 206512 && (lumi >= 203 && lumi <= 212))
    return true;
  if (run == 206512 && (lumi >= 214 && lumi <= 332))
    return true;
  if (run == 206512 && (lumi >= 334 && lumi <= 584))
    return true;
  if (run == 206512 && (lumi >= 587 && lumi <= 604))
    return true;
  if (run == 206512 && (lumi >= 607 && lumi <= 1005))
    return true;
  if (run == 206512 && (lumi >= 1008 && lumi <= 1123))
    return true;
  if (run == 206512 && (lumi >= 1126 && lumi <= 1163))
    return true;
  if (run == 206512 && (lumi >= 1165 && lumi <= 1211))
    return true;
  if (run == 206513 && (lumi >= 3 && lumi <= 39))
    return true;
  if (run == 206513 && (lumi >= 42 && lumi <= 188))
    return true;
  if (run == 206513 && (lumi >= 191 && lumi <= 234))
    return true;
  if (run == 206513 && (lumi >= 237 && lumi <= 238))
    return true;
  if (run == 206513 && (lumi >= 241 && lumi <= 323))
    return true;
  if (run == 207454 && (lumi >= 79 && lumi <= 95))
    return true;
  if (run == 207454 && (lumi >= 98 && lumi <= 123))
    return true;
  if (run == 207454 && (lumi >= 126 && lumi <= 259))
    return true;
  if (run == 207454 && (lumi >= 261 && lumi <= 363))
    return true;
  if (run == 207454 && (lumi >= 365 && lumi <= 458))
    return true;
  if (run == 207454 && (lumi >= 461 && lumi <= 498))
    return true;
  if (run == 207454 && (lumi >= 501 && lumi <= 609))
    return true;
  if (run == 207454 && (lumi >= 612 && lumi <= 632))
    return true;
  if (run == 207454 && (lumi >= 635 && lumi <= 781))
    return true;
  if (run == 207454 && (lumi >= 784 && lumi <= 866))
    return true;
  if (run == 207454 && (lumi >= 869 && lumi <= 974))
    return true;
  if (run == 207454 && (lumi >= 977 && lumi <= 1064))
    return true;
  if (run == 207454 && (lumi >= 1067 && lumi <= 1079))
    return true;
  if (run == 207454 && (lumi >= 1081 && lumi <= 1321))
    return true;
  if (run == 207454 && (lumi >= 1323 && lumi <= 1464))
    return true;
  if (run == 207454 && (lumi >= 1467 && lumi <= 1569))
    return true;
  if (run == 207454 && (lumi >= 1571 && lumi <= 1604))
    return true;
  if (run == 207454 && (lumi >= 1607 && lumi <= 1712))
    return true;
  if (run == 207454 && (lumi >= 1714 && lumi <= 1988))
    return true;
  if (run == 198372 && (lumi >= 57 && lumi <= 110))
    return true;
  if (run == 194223 && (lumi >= 61 && lumi <= 112))
    return true;
  if (run == 200600 && (lumi >= 1 && lumi <= 183))
    return true;
  if (run == 200600 && (lumi >= 186 && lumi <= 299))
    return true;
  if (run == 200600 && (lumi >= 302 && lumi <= 313))
    return true;
  if (run == 200600 && (lumi >= 316 && lumi <= 324))
    return true;
  if (run == 200600 && (lumi >= 327 && lumi <= 334))
    return true;
  if (run == 200600 && (lumi >= 336 && lumi <= 397))
    return true;
  if (run == 200600 && (lumi >= 399 && lumi <= 417))
    return true;
  if (run == 200600 && (lumi >= 420 && lumi <= 526))
    return true;
  if (run == 200600 && (lumi >= 529 && lumi <= 591))
    return true;
  if (run == 200600 && (lumi >= 594 && lumi <= 596))
    return true;
  if (run == 200600 && (lumi >= 598 && lumi <= 609))
    return true;
  if (run == 200600 && (lumi >= 611 && lumi <= 660))
    return true;
  if (run == 200600 && (lumi >= 663 && lumi <= 823))
    return true;
  if (run == 200600 && (lumi >= 826 && lumi <= 900))
    return true;
  if (run == 200600 && (lumi >= 902 && lumi <= 943))
    return true;
  if (run == 200600 && (lumi >= 945 && lumi <= 1139))
    return true;
  if (run == 195749 && (lumi >= 1 && lumi <= 8))
    return true;
  if (run == 195749 && (lumi >= 10 && lumi <= 33))
    return true;
  if (run == 195749 && (lumi >= 36 && lumi <= 131))
    return true;
  if (run == 203830 && (lumi >= 82 && lumi <= 182))
    return true;
  if (run == 201173 && (lumi >= 1 && lumi <= 194))
    return true;
  if (run == 201173 && (lumi >= 197 && lumi <= 586))
    return true;
  if (run == 201174 && (lumi >= 1 && lumi <= 214))
    return true;
  if (run == 201174 && (lumi >= 216 && lumi <= 263))
    return true;
  if (run == 201174 && (lumi >= 265 && lumi <= 339))
    return true;
  if (run == 201174 && (lumi >= 342 && lumi <= 451))
    return true;
  if (run == 195013 && (lumi >= 94 && lumi <= 144))
    return true;
  if (run == 195013 && (lumi >= 146 && lumi <= 185))
    return true;
  if (run == 195013 && (lumi >= 187 && lumi <= 206))
    return true;
  if (run == 195013 && (lumi >= 208 && lumi <= 299))
    return true;
  if (run == 195013 && (lumi >= 302 && lumi <= 324))
    return true;
  if (run == 195013 && (lumi >= 326 && lumi <= 366))
    return true;
  if (run == 195013 && (lumi >= 369 && lumi <= 447))
    return true;
  if (run == 195013 && (lumi >= 450 && lumi <= 526))
    return true;
  if (run == 195013 && (lumi >= 528 && lumi <= 541))
    return true;
  if (run == 195014 && (lumi >= 1 && lumi <= 6))
    return true;
  if (run == 195014 && (lumi >= 9 && lumi <= 119))
    return true;
  if (run == 195014 && (lumi >= 121 && lumi <= 148))
    return true;
  if (run == 195015 && (lumi >= 1 && lumi <= 13))
    return true;
  if (run == 195016 && (lumi >= 1 && lumi <= 21))
    return true;
  if (run == 195016 && (lumi >= 23 && lumi <= 55))
    return true;
  if (run == 195016 && (lumi >= 58 && lumi <= 63))
    return true;
  if (run == 195016 && (lumi >= 65 && lumi <= 174))
    return true;
  if (run == 195016 && (lumi >= 177 && lumi <= 184))
    return true;
  if (run == 195016 && (lumi >= 186 && lumi <= 241))
    return true;
  if (run == 195016 && (lumi >= 243 && lumi <= 246))
    return true;
  if (run == 195016 && (lumi >= 248 && lumi <= 251))
    return true;
  if (run == 195016 && (lumi >= 254 && lumi <= 367))
    return true;
  if (run == 195016 && (lumi >= 370 && lumi <= 422))
    return true;
  if (run == 195016 && (lumi >= 425 && lumi <= 560))
    return true;
  if (run == 195016 && (lumi >= 563 && lumi <= 569))
    return true;
  if (run == 198207 && (lumi >= 1 && lumi <= 97))
    return true;
  if (run == 201197 && (lumi >= 1 && lumi <= 23))
    return true;
  if (run == 201196 && (lumi >= 1 && lumi <= 238))
    return true;
  if (run == 201196 && (lumi >= 241 && lumi <= 278))
    return true;
  if (run == 201196 && (lumi >= 286 && lumi <= 299))
    return true;
  if (run == 201196 && (lumi >= 302 && lumi <= 338))
    return true;
  if (run == 201196 && (lumi >= 341 && lumi <= 515))
    return true;
  if (run == 201196 && (lumi >= 518 && lumi <= 720))
    return true;
  if (run == 201196 && (lumi >= 723 && lumi <= 789))
    return true;
  if (run == 201196 && (lumi >= 803 && lumi <= 841))
    return true;
  if (run == 196495 && (lumi >= 114 && lumi <= 180))
    return true;
  if (run == 196495 && (lumi >= 182 && lumi <= 272))
    return true;
  if (run == 191248 && (lumi >= 1 && lumi <= 103))
    return true;
  if (run == 201193 && (lumi >= 1 && lumi <= 19))
    return true;
  if (run == 191247 && (lumi >= 1 && lumi <= 153))
    return true;
  if (run == 191247 && (lumi >= 156 && lumi <= 280))
    return true;
  if (run == 191247 && (lumi >= 283 && lumi <= 606))
    return true;
  if (run == 191247 && (lumi >= 608 && lumi <= 620))
    return true;
  if (run == 191247 && (lumi >= 622 && lumi <= 818))
    return true;
  if (run == 191247 && (lumi >= 821 && lumi <= 834))
    return true;
  if (run == 191247 && (lumi >= 837 && lumi <= 1031))
    return true;
  if (run == 191247 && (lumi >= 1034 && lumi <= 1046))
    return true;
  if (run == 191247 && (lumi >= 1049 && lumi <= 1140))
    return true;
  if (run == 191247 && (lumi >= 1143 && lumi <= 1187))
    return true;
  if (run == 191247 && (lumi >= 1190 && lumi <= 1214))
    return true;
  if (run == 191247 && (lumi >= 1217 && lumi <= 1224))
    return true;
  if (run == 198208 && (lumi >= 1 && lumi <= 92))
    return true;
  if (run == 198208 && (lumi >= 94 && lumi <= 134))
    return true;
  if (run == 198208 && (lumi >= 137 && lumi <= 147))
    return true;
  if (run == 198208 && (lumi >= 150 && lumi <= 209))
    return true;
  if (run == 202237 && (lumi >= 39 && lumi <= 128))
    return true;
  if (run == 202237 && (lumi >= 131 && lumi <= 131))
    return true;
  if (run == 202237 && (lumi >= 134 && lumi <= 219))
    return true;
  if (run == 202237 && (lumi >= 222 && lumi <= 235))
    return true;
  if (run == 202237 && (lumi >= 238 && lumi <= 275))
    return true;
  if (run == 202237 && (lumi >= 277 && lumi <= 289))
    return true;
  if (run == 202237 && (lumi >= 291 && lumi <= 316))
    return true;
  if (run == 202237 && (lumi >= 319 && lumi <= 419))
    return true;
  if (run == 202237 && (lumi >= 422 && lumi <= 538))
    return true;
  if (run == 202237 && (lumi >= 540 && lumi <= 936))
    return true;
  if (run == 202237 && (lumi >= 939 && lumi <= 950))
    return true;
  if (run == 202237 && (lumi >= 952 && lumi <= 976))
    return true;
  if (run == 202237 && (lumi >= 979 && lumi <= 1079))
    return true;
  if (run == 199008 && (lumi >= 75 && lumi <= 93))
    return true;
  if (run == 199008 && (lumi >= 95 && lumi <= 121))
    return true;
  if (run == 199008 && (lumi >= 124 && lumi <= 208))
    return true;
  if (run == 199008 && (lumi >= 211 && lumi <= 331))
    return true;
  if (run == 199008 && (lumi >= 333 && lumi <= 373))
    return true;
  if (run == 199008 && (lumi >= 376 && lumi <= 482))
    return true;
  if (run == 199008 && (lumi >= 485 && lumi <= 605))
    return true;
  if (run == 199008 && (lumi >= 608 && lumi <= 644))
    return true;
  if (run == 198272 && (lumi >= 1 && lumi <= 185))
    return true;
  if (run == 198272 && (lumi >= 188 && lumi <= 245))
    return true;
  if (run == 198272 && (lumi >= 248 && lumi <= 314))
    return true;
  if (run == 198272 && (lumi >= 317 && lumi <= 433))
    return true;
  if (run == 198272 && (lumi >= 436 && lumi <= 444))
    return true;
  if (run == 198272 && (lumi >= 454 && lumi <= 620))
    return true;
  if (run == 198271 && (lumi >= 1 && lumi <= 91))
    return true;
  if (run == 198271 && (lumi >= 93 && lumi <= 170))
    return true;
  if (run == 198271 && (lumi >= 173 && lumi <= 299))
    return true;
  if (run == 198271 && (lumi >= 301 && lumi <= 450))
    return true;
  if (run == 198271 && (lumi >= 453 && lumi <= 513))
    return true;
  if (run == 198271 && (lumi >= 516 && lumi <= 616))
    return true;
  if (run == 198271 && (lumi >= 619 && lumi <= 628))
    return true;
  if (run == 198271 && (lumi >= 631 && lumi <= 791))
    return true;
  if (run == 198271 && (lumi >= 793 && lumi <= 797))
    return true;
  if (run == 205694 && (lumi >= 1 && lumi <= 205))
    return true;
  if (run == 205694 && (lumi >= 208 && lumi <= 230))
    return true;
  if (run == 205694 && (lumi >= 233 && lumi <= 347))
    return true;
  if (run == 205694 && (lumi >= 350 && lumi <= 452))
    return true;
  if (run == 205694 && (lumi >= 455 && lumi <= 593))
    return true;
  if (run == 205694 && (lumi >= 595 && lumi <= 890))
    return true;
  if (run == 205690 && (lumi >= 1 && lumi <= 40))
    return true;
  if (run == 195841 && (lumi >= 74 && lumi <= 85))
    return true;
  if (run == 206258 && (lumi >= 1 && lumi <= 36))
    return true;
  if (run == 206258 && (lumi >= 39 && lumi <= 223))
    return true;
  if (run == 206258 && (lumi >= 226 && lumi <= 249))
    return true;
  if (run == 206257 && (lumi >= 1 && lumi <= 29))
    return true;
  if (run == 202044 && (lumi >= 84 && lumi <= 101))
    return true;
  if (run == 202044 && (lumi >= 104 && lumi <= 266))
    return true;
  if (run == 202044 && (lumi >= 268 && lumi <= 461))
    return true;
  if (run == 202044 && (lumi >= 463 && lumi <= 466))
    return true;
  if (run == 202045 && (lumi >= 1 && lumi <= 30))
    return true;
  if (run == 202045 && (lumi >= 33 && lumi <= 72))
    return true;
  if (run == 202045 && (lumi >= 75 && lumi <= 528))
    return true;
  if (run == 202045 && (lumi >= 531 && lumi <= 601))
    return true;
  if (run == 202045 && (lumi >= 603 && lumi <= 785))
    return true;
  if (run == 202045 && (lumi >= 788 && lumi <= 809))
    return true;
  if (run == 202045 && (lumi >= 822 && lumi <= 823))
    return true;
  if (run == 200180 && (lumi >= 1 && lumi <= 18))
    return true;
  if (run == 200186 && (lumi >= 1 && lumi <= 3))
    return true;
  if (run == 200186 && (lumi >= 6 && lumi <= 24))
    return true;
  if (run == 194711 && (lumi >= 1 && lumi <= 7))
    return true;
  if (run == 194711 && (lumi >= 9 && lumi <= 619))
    return true;
  if (run == 194712 && (lumi >= 1 && lumi <= 56))
    return true;
  if (run == 194712 && (lumi >= 61 && lumi <= 418))
    return true;
  if (run == 194712 && (lumi >= 420 && lumi <= 625))
    return true;
  if (run == 194712 && (lumi >= 627 && lumi <= 759))
    return true;
  if (run == 205826 && (lumi >= 80 && lumi <= 232))
    return true;
  if (run == 205826 && (lumi >= 235 && lumi <= 303))
    return true;
  if (run == 205826 && (lumi >= 306 && lumi <= 468))
    return true;
  if (run == 191367 && (lumi >= 1 && lumi <= 2))
    return true;
  if (run == 206744 && (lumi >= 49 && lumi <= 157))
    return true;
  if (run == 206744 && (lumi >= 160 && lumi <= 192))
    return true;
  if (run == 206744 && (lumi >= 195 && lumi <= 395))
    return true;
  if (run == 206744 && (lumi >= 398 && lumi <= 452))
    return true;
  if (run == 201202 && (lumi >= 1 && lumi <= 437))
    return true;
  if (run == 205217 && (lumi >= 1 && lumi <= 12))
    return true;
  if (run == 205217 && (lumi >= 16 && lumi <= 111))
    return true;
  if (run == 205217 && (lumi >= 113 && lumi <= 171))
    return true;
  if (run == 205217 && (lumi >= 174 && lumi <= 250))
    return true;
  if (run == 205217 && (lumi >= 253 && lumi <= 318))
    return true;
  if (run == 195937 && (lumi >= 1 && lumi <= 28))
    return true;
  if (run == 195937 && (lumi >= 31 && lumi <= 186))
    return true;
  if (run == 195937 && (lumi >= 188 && lumi <= 396))
    return true;
  if (run == 195304 && (lumi >= 1 && lumi <= 3))
    return true;
  if (run == 195304 && (lumi >= 6 && lumi <= 22))
    return true;
  if (run == 195304 && (lumi >= 27 && lumi <= 80))
    return true;
  if (run == 195304 && (lumi >= 83 && lumi <= 100))
    return true;
  if (run == 195304 && (lumi >= 103 && lumi <= 154))
    return true;
  if (run == 195304 && (lumi >= 157 && lumi <= 341))
    return true;
  if (run == 195304 && (lumi >= 344 && lumi <= 588))
    return true;
  if (run == 195304 && (lumi >= 590 && lumi <= 727))
    return true;
  if (run == 195304 && (lumi >= 729 && lumi <= 1003))
    return true;
  if (run == 195304 && (lumi >= 1006 && lumi <= 1079))
    return true;
  if (run == 195304 && (lumi >= 1083 && lumi <= 1140))
    return true;
  if (run == 195304 && (lumi >= 1143 && lumi <= 1229))
    return true;
  if (run == 205526 && (lumi >= 1 && lumi <= 269))
    return true;
  if (run == 205526 && (lumi >= 272 && lumi <= 277))
    return true;
  if (run == 205526 && (lumi >= 280 && lumi <= 332))
    return true;
  if (run == 195930 && (lumi >= 1 && lumi <= 77))
    return true;
  if (run == 195930 && (lumi >= 80 && lumi <= 176))
    return true;
  if (run == 195930 && (lumi >= 179 && lumi <= 526))
    return true;
  if (run == 195930 && (lumi >= 529 && lumi <= 596))
    return true;
  if (run == 199754 && (lumi >= 1 && lumi <= 203))
    return true;
  if (run == 199754 && (lumi >= 205 && lumi <= 325))
    return true;
  if (run == 199754 && (lumi >= 328 && lumi <= 457))
    return true;
  if (run == 199754 && (lumi >= 459 && lumi <= 607))
    return true;
  if (run == 199754 && (lumi >= 610 && lumi <= 613))
    return true;
  if (run == 199754 && (lumi >= 615 && lumi <= 806))
    return true;
  if (run == 199754 && (lumi >= 808 && lumi <= 998))
    return true;
  if (run == 199753 && (lumi >= 1 && lumi <= 59))
    return true;
  if (run == 199752 && (lumi >= 1 && lumi <= 141))
    return true;
  if (run == 199752 && (lumi >= 144 && lumi <= 180))
    return true;
  if (run == 199752 && (lumi >= 182 && lumi <= 186))
    return true;
  if (run == 199752 && (lumi >= 188 && lumi <= 211))
    return true;
  if (run == 199752 && (lumi >= 214 && lumi <= 322))
    return true;
  if (run == 199698 && (lumi >= 72 && lumi <= 94))
    return true;
  if (run == 199698 && (lumi >= 96 && lumi <= 127))
    return true;
  if (run == 199699 && (lumi >= 1 && lumi <= 154))
    return true;
  if (run == 199699 && (lumi >= 157 && lumi <= 169))
    return true;
  if (run == 199699 && (lumi >= 172 && lumi <= 410))
    return true;
  if (run == 199699 && (lumi >= 412 && lumi <= 756))
    return true;
  if (run == 191062 && (lumi >= 1 && lumi <= 1))
    return true;
  if (run == 191062 && (lumi >= 3 && lumi <= 3))
    return true;
  if (run == 191062 && (lumi >= 5 && lumi <= 214))
    return true;
  if (run == 191062 && (lumi >= 216 && lumi <= 541))
    return true;
  if (run == 193207 && (lumi >= 54 && lumi <= 182))
    return true;
  if (run == 190895 && (lumi >= 64 && lumi <= 202))
    return true;
  if (run == 190895 && (lumi >= 210 && lumi <= 302))
    return true;
  if (run == 190895 && (lumi >= 305 && lumi <= 584))
    return true;
  if (run == 190895 && (lumi >= 587 && lumi <= 948))
    return true;
  if (run == 194199 && (lumi >= 96 && lumi <= 227))
    return true;
  if (run == 194199 && (lumi >= 229 && lumi <= 336))
    return true;
  if (run == 194199 && (lumi >= 339 && lumi <= 402))
    return true;
  if (run == 196334 && (lumi >= 59 && lumi <= 111))
    return true;
  if (run == 196334 && (lumi >= 113 && lumi <= 123))
    return true;
  if (run == 196334 && (lumi >= 126 && lumi <= 132))
    return true;
  if (run == 196334 && (lumi >= 135 && lumi <= 167))
    return true;
  if (run == 196334 && (lumi >= 170 && lumi <= 193))
    return true;
  if (run == 196334 && (lumi >= 196 && lumi <= 257))
    return true;
  if (run == 196334 && (lumi >= 259 && lumi <= 267))
    return true;
  if (run == 196334 && (lumi >= 270 && lumi <= 289))
    return true;
  if (run == 196334 && (lumi >= 292 && lumi <= 342))
    return true;
  if (run == 199356 && (lumi >= 95 && lumi <= 121))
    return true;
  if (run == 199356 && (lumi >= 123 && lumi <= 168))
    return true;
  if (run == 199356 && (lumi >= 171 && lumi <= 205))
    return true;
  if (run == 199356 && (lumi >= 208 && lumi <= 231))
    return true;
  if (run == 195929 && (lumi >= 1 && lumi <= 29))
    return true;
  if (run == 195303 && (lumi >= 109 && lumi <= 191))
    return true;
  if (run == 195303 && (lumi >= 194 && lumi <= 277))
    return true;
  if (run == 195303 && (lumi >= 280 && lumi <= 310))
    return true;
  if (run == 195303 && (lumi >= 312 && lumi <= 316))
    return true;
  if (run == 195303 && (lumi >= 318 && lumi <= 409))
    return true;
  if (run == 206550 && (lumi >= 77 && lumi <= 132))
    return true;
  if (run == 206550 && (lumi >= 135 && lumi <= 144))
    return true;
  if (run == 207490 && (lumi >= 1 && lumi <= 48))
    return true;
  if (run == 207490 && (lumi >= 51 && lumi <= 111))
    return true;
  if (run == 207491 && (lumi >= 1 && lumi <= 176))
    return true;
  if (run == 207491 && (lumi >= 179 && lumi <= 458))
    return true;
  if (run == 206389 && (lumi >= 88 && lumi <= 185))
    return true;
  if (run == 206389 && (lumi >= 187 && lumi <= 249))
    return true;
  if (run == 206389 && (lumi >= 252 && lumi <= 272))
    return true;
  if (run == 206389 && (lumi >= 275 && lumi <= 392))
    return true;
  if (run == 203987 && (lumi >= 1 && lumi <= 9))
    return true;
  if (run == 203987 && (lumi >= 12 && lumi <= 241))
    return true;
  if (run == 203987 && (lumi >= 243 && lumi <= 339))
    return true;
  if (run == 203987 && (lumi >= 342 && lumi <= 781))
    return true;
  if (run == 203987 && (lumi >= 784 && lumi <= 1014))
    return true;
  if (run == 194789 && (lumi >= 1 && lumi <= 18))
    return true;
  if (run == 194789 && (lumi >= 21 && lumi <= 32))
    return true;
  if (run == 194789 && (lumi >= 34 && lumi <= 80))
    return true;
  if (run == 194789 && (lumi >= 82 && lumi <= 166))
    return true;
  if (run == 194789 && (lumi >= 168 && lumi <= 269))
    return true;
  if (run == 194789 && (lumi >= 272 && lumi <= 405))
    return true;
  if (run == 194789 && (lumi >= 409 && lumi <= 414))
    return true;
  if (run == 194789 && (lumi >= 417 && lumi <= 427))
    return true;
  if (run == 194789 && (lumi >= 430 && lumi <= 566))
    return true;
  if (run == 208540 && (lumi >= 1 && lumi <= 26))
    return true;
  if (run == 208540 && (lumi >= 29 && lumi <= 98))
    return true;
  if (run == 204101 && (lumi >= 1 && lumi <= 74))
    return true;
  if (run == 204100 && (lumi >= 117 && lumi <= 139))
    return true;
  if (run == 202469 && (lumi >= 87 && lumi <= 158))
    return true;
  if (run == 202469 && (lumi >= 160 && lumi <= 174))
    return true;
  if (run == 202469 && (lumi >= 177 && lumi <= 352))
    return true;
  if (run == 207924 && (lumi >= 1 && lumi <= 52))
    return true;
  if (run == 207924 && (lumi >= 54 && lumi <= 171))
    return true;
  if (run == 207924 && (lumi >= 173 && lumi <= 178))
    return true;
  if (run == 207924 && (lumi >= 181 && lumi <= 339))
    return true;
  if (run == 207922 && (lumi >= 1 && lumi <= 69))
    return true;
  if (run == 207922 && (lumi >= 71 && lumi <= 100))
    return true;
  if (run == 207922 && (lumi >= 103 && lumi <= 126))
    return true;
  if (run == 207922 && (lumi >= 129 && lumi <= 242))
    return true;
  if (run == 207922 && (lumi >= 274 && lumi <= 291))
    return true;
  if (run == 207920 && (lumi >= 84 && lumi <= 146))
    return true;
  if (run == 207920 && (lumi >= 149 && lumi <= 241))
    return true;
  if (run == 207920 && (lumi >= 243 && lumi <= 261))
    return true;
  if (run == 207920 && (lumi >= 264 && lumi <= 291))
    return true;
  if (run == 207920 && (lumi >= 294 && lumi <= 486))
    return true;
  if (run == 207920 && (lumi >= 489 && lumi <= 518))
    return true;
  if (run == 207920 && (lumi >= 520 && lumi <= 598))
    return true;
  if (run == 207920 && (lumi >= 600 && lumi <= 708))
    return true;
  if (run == 207920 && (lumi >= 710 && lumi <= 826))
    return true;
  if (run == 207921 && (lumi >= 1 && lumi <= 37))
    return true;
  if (run == 207921 && (lumi >= 40 && lumi <= 58))
    return true;
  if (run == 194224 && (lumi >= 1 && lumi <= 126))
    return true;
  if (run == 194224 && (lumi >= 129 && lumi <= 206))
    return true;
  if (run == 194224 && (lumi >= 208 && lumi <= 250))
    return true;
  if (run == 194224 && (lumi >= 253 && lumi <= 309))
    return true;
  if (run == 194224 && (lumi >= 312 && lumi <= 386))
    return true;
  if (run == 194224 && (lumi >= 389 && lumi <= 412))
    return true;
  if (run == 194225 && (lumi >= 1 && lumi <= 23))
    return true;
  if (run == 194225 && (lumi >= 26 && lumi <= 47))
    return true;
  if (run == 194225 && (lumi >= 49 && lumi <= 85))
    return true;
  if (run == 194225 && (lumi >= 88 && lumi <= 149))
    return true;
  if (run == 196252 && (lumi >= 1 && lumi <= 35))
    return true;
  if (run == 196250 && (lumi >= 1 && lumi <= 2))
    return true;
  if (run == 196250 && (lumi >= 5 && lumi <= 265))
    return true;
  if (run == 196250 && (lumi >= 267 && lumi <= 426))
    return true;
  if (run == 204564 && (lumi >= 1 && lumi <= 84))
    return true;
  if (run == 204564 && (lumi >= 87 && lumi <= 89))
    return true;
  if (run == 204564 && (lumi >= 92 && lumi <= 159))
    return true;
  if (run == 204564 && (lumi >= 161 && lumi <= 187))
    return true;
  if (run == 204564 && (lumi >= 190 && lumi <= 191))
    return true;
  if (run == 204564 && (lumi >= 193 && lumi <= 293))
    return true;
  if (run == 204564 && (lumi >= 296 && lumi <= 315))
    return true;
  if (run == 204564 && (lumi >= 317 && lumi <= 340))
    return true;
  if (run == 204564 && (lumi >= 343 && lumi <= 427))
    return true;
  if (run == 204564 && (lumi >= 429 && lumi <= 434))
    return true;
  if (run == 204564 && (lumi >= 437 && lumi <= 735))
    return true;
  if (run == 204564 && (lumi >= 737 && lumi <= 855))
    return true;
  if (run == 204564 && (lumi >= 858 && lumi <= 1206))
    return true;
  if (run == 204564 && (lumi >= 1209 && lumi <= 1248))
    return true;
  if (run == 204564 && (lumi >= 1251 && lumi <= 1284))
    return true;
  if (run == 201708 && (lumi >= 1 && lumi <= 79))
    return true;
  if (run == 208538 && (lumi >= 2 && lumi <= 43))
    return true;
  if (run == 201705 && (lumi >= 65 && lumi <= 73))
    return true;
  if (run == 201705 && (lumi >= 75 && lumi <= 109))
    return true;
  if (run == 201705 && (lumi >= 111 && lumi <= 187))
    return true;
  if (run == 201706 && (lumi >= 1 && lumi <= 62))
    return true;
  if (run == 201707 && (lumi >= 1 && lumi <= 23))
    return true;
  if (run == 201707 && (lumi >= 26 && lumi <= 42))
    return true;
  if (run == 201707 && (lumi >= 45 && lumi <= 115))
    return true;
  if (run == 201707 && (lumi >= 118 && lumi <= 130))
    return true;
  if (run == 201707 && (lumi >= 133 && lumi <= 160))
    return true;
  if (run == 201707 && (lumi >= 163 && lumi <= 276))
    return true;
  if (run == 201707 && (lumi >= 279 && lumi <= 471))
    return true;
  if (run == 201707 && (lumi >= 473 && lumi <= 511))
    return true;
  if (run == 201707 && (lumi >= 514 && lumi <= 545))
    return true;
  if (run == 201707 && (lumi >= 547 && lumi <= 570))
    return true;
  if (run == 201707 && (lumi >= 572 && lumi <= 622))
    return true;
  if (run == 201707 && (lumi >= 625 && lumi <= 735))
    return true;
  if (run == 201707 && (lumi >= 738 && lumi <= 806))
    return true;
  if (run == 201707 && (lumi >= 809 && lumi <= 876))
    return true;
  if (run == 201707 && (lumi >= 879 && lumi <= 964))
    return true;
  if (run == 198249 && (lumi >= 1 && lumi <= 7))
    return true;
  if (run == 204576 && (lumi >= 49 && lumi <= 192))
    return true;
  if (run == 204576 && (lumi >= 195 && lumi <= 301))
    return true;
  if (run == 201097 && (lumi >= 83 && lumi <= 136))
    return true;
  if (run == 201097 && (lumi >= 138 && lumi <= 245))
    return true;
  if (run == 201097 && (lumi >= 248 && lumi <= 300))
    return true;
  if (run == 201097 && (lumi >= 303 && lumi <= 370))
    return true;
  if (run == 201097 && (lumi >= 372 && lumi <= 429))
    return true;
  if (run == 201097 && (lumi >= 432 && lumi <= 497))
    return true;
  if (run == 199318 && (lumi >= 65 && lumi <= 138))
    return true;
  if (run == 195656 && (lumi >= 1 && lumi <= 362))
    return true;
  if (run == 194533 && (lumi >= 46 && lumi <= 379))
    return true;
  if (run == 194533 && (lumi >= 382 && lumi <= 415))
    return true;
  if (run == 194533 && (lumi >= 417 && lumi <= 618))
    return true;
  if (run == 194533 && (lumi >= 620 && lumi <= 872))
    return true;
  if (run == 195658 && (lumi >= 1 && lumi <= 37))
    return true;
  if (run == 195658 && (lumi >= 40 && lumi <= 362))
    return true;
  if (run == 195658 && (lumi >= 364 && lumi <= 382))
    return true;
  if (run == 195658 && (lumi >= 384 && lumi <= 386))
    return true;
  if (run == 196027 && (lumi >= 1 && lumi <= 55))
    return true;
  if (run == 196027 && (lumi >= 58 && lumi <= 119))
    return true;
  if (run == 196027 && (lumi >= 121 && lumi <= 155))
    return true;
  if (run == 196027 && (lumi >= 158 && lumi <= 186))
    return true;
  if (run == 201794 && (lumi >= 58 && lumi <= 94))
    return true;
  if (run == 200381 && (lumi >= 8 && lumi <= 15))
    return true;
  if (run == 200381 && (lumi >= 18 && lumi <= 36))
    return true;
  if (run == 200381 && (lumi >= 38 && lumi <= 89))
    return true;
  if (run == 200381 && (lumi >= 91 && lumi <= 195))
    return true;
  if (run == 191721 && (lumi >= 1 && lumi <= 1))
    return true;
  if (run == 191721 && (lumi >= 3 && lumi <= 34))
    return true;
  if (run == 191721 && (lumi >= 36 && lumi <= 183))
    return true;
  if (run == 191721 && (lumi >= 186 && lumi <= 189))
    return true;
  if (run == 191726 && (lumi >= 1 && lumi <= 13))
    return true;
  if (run == 202088 && (lumi >= 1 && lumi <= 286))
    return true;
  if (run == 202087 && (lumi >= 1 && lumi <= 25))
    return true;
  if (run == 202087 && (lumi >= 28 && lumi <= 208))
    return true;
  if (run == 202087 && (lumi >= 210 && lumi <= 357))
    return true;
  if (run == 202087 && (lumi >= 359 && lumi <= 652))
    return true;
  if (run == 202087 && (lumi >= 655 && lumi <= 853))
    return true;
  if (run == 202087 && (lumi >= 856 && lumi <= 1093))
    return true;
  if (run == 202084 && (lumi >= 83 && lumi <= 156))
    return true;
  if (run == 202084 && (lumi >= 159 && lumi <= 177))
    return true;
  if (run == 202084 && (lumi >= 179 && lumi <= 180))
    return true;
  if (run == 202084 && (lumi >= 182 && lumi <= 239))
    return true;
  if (run == 194120 && (lumi >= 1 && lumi <= 162))
    return true;
  if (run == 194120 && (lumi >= 165 && lumi <= 406))
    return true;
  if (run == 208390 && (lumi >= 72 && lumi <= 128))
    return true;
  if (run == 208390 && (lumi >= 130 && lumi <= 169))
    return true;
  if (run == 208391 && (lumi >= 52 && lumi <= 82))
    return true;
  if (run == 208391 && (lumi >= 84 && lumi <= 162))
    return true;
  if (run == 208391 && (lumi >= 164 && lumi <= 216))
    return true;
  if (run == 208391 && (lumi >= 219 && lumi <= 493))
    return true;
  if (run == 208391 && (lumi >= 495 && lumi <= 498))
    return true;
  if (run == 208391 && (lumi >= 500 && lumi <= 523))
    return true;
  if (run == 208391 && (lumi >= 526 && lumi <= 533))
    return true;
  if (run == 208391 && (lumi >= 535 && lumi <= 588))
    return true;
  if (run == 208391 && (lumi >= 591 && lumi <= 660))
    return true;
  if (run == 208391 && (lumi >= 663 && lumi <= 869))
    return true;
  if (run == 194455 && (lumi >= 45 && lumi <= 64))
    return true;
  if (run == 194455 && (lumi >= 67 && lumi <= 140))
    return true;
  if (run == 194455 && (lumi >= 142 && lumi <= 255))
    return true;
  if (run == 194455 && (lumi >= 293 && lumi <= 303))
    return true;
  if (run == 196509 && (lumi >= 1 && lumi <= 68))
    return true;
  if (run == 195530 && (lumi >= 1 && lumi <= 80))
    return true;
  if (run == 195530 && (lumi >= 82 && lumi <= 104))
    return true;
  if (run == 195530 && (lumi >= 107 && lumi <= 156))
    return true;
  if (run == 195530 && (lumi >= 159 && lumi <= 300))
    return true;
  if (run == 195530 && (lumi >= 302 && lumi <= 405))
    return true;
  if (run == 205620 && (lumi >= 1 && lumi <= 175))
    return true;
  if (run == 207100 && (lumi >= 1 && lumi <= 91))
    return true;
  if (run == 207100 && (lumi >= 94 && lumi <= 94))
    return true;
  if (run == 205781 && (lumi >= 1 && lumi <= 89))
    return true;
  if (run == 205781 && (lumi >= 91 && lumi <= 197))
    return true;
  if (run == 205781 && (lumi >= 200 && lumi <= 502))
    return true;
  if (run == 191201 && (lumi >= 38 && lumi <= 49))
    return true;
  if (run == 191201 && (lumi >= 52 && lumi <= 79))
    return true;
  if (run == 202305 && (lumi >= 1 && lumi <= 89))
    return true;
  if (run == 202305 && (lumi >= 92 && lumi <= 130))
    return true;
  if (run == 202305 && (lumi >= 133 && lumi <= 323))
    return true;
  if (run == 191202 && (lumi >= 1 && lumi <= 64))
    return true;
  if (run == 191202 && (lumi >= 66 && lumi <= 68))
    return true;
  if (run == 191202 && (lumi >= 87 && lumi <= 105))
    return true;
  if (run == 191202 && (lumi >= 108 && lumi <= 118))
    return true;
  if (run == 193557 && (lumi >= 1 && lumi <= 84))
    return true;
  if (run == 193556 && (lumi >= 41 && lumi <= 83))
    return true;
  if (run == 195970 && (lumi >= 44 && lumi <= 49))
    return true;
  if (run == 195970 && (lumi >= 51 && lumi <= 85))
    return true;
  if (run == 207882 && (lumi >= 22 && lumi <= 45))
    return true;
  if (run == 200161 && (lumi >= 1 && lumi <= 97))
    return true;
  if (run == 200161 && (lumi >= 100 && lumi <= 112))
    return true;
  if (run == 200160 && (lumi >= 52 && lumi <= 68))
    return true;
  if (run == 202014 && (lumi >= 1 && lumi <= 5))
    return true;
  if (run == 202014 && (lumi >= 8 && lumi <= 14))
    return true;
  if (run == 202014 && (lumi >= 16 && lumi <= 18))
    return true;
  if (run == 202014 && (lumi >= 20 && lumi <= 77))
    return true;
  if (run == 202014 && (lumi >= 79 && lumi <= 102))
    return true;
  if (run == 202014 && (lumi >= 104 && lumi <= 174))
    return true;
  if (run == 202014 && (lumi >= 177 && lumi <= 190))
    return true;
  if (run == 202014 && (lumi >= 192 && lumi <= 196))
    return true;
  if (run == 195655 && (lumi >= 1 && lumi <= 129))
    return true;
  if (run == 195655 && (lumi >= 131 && lumi <= 184))
    return true;
  if (run == 195655 && (lumi >= 186 && lumi <= 260))
    return true;
  if (run == 195655 && (lumi >= 263 && lumi <= 350))
    return true;
  if (run == 195655 && (lumi >= 353 && lumi <= 446))
    return true;
  if (run == 195655 && (lumi >= 448 && lumi <= 483))
    return true;
  if (run == 195655 && (lumi >= 485 && lumi <= 498))
    return true;
  if (run == 208307 && (lumi >= 2 && lumi <= 42))
    return true;
  if (run == 208307 && (lumi >= 45 && lumi <= 45))
    return true;
  if (run == 208307 && (lumi >= 47 && lumi <= 70))
    return true;
  if (run == 208307 && (lumi >= 72 && lumi <= 147))
    return true;
  if (run == 208307 && (lumi >= 150 && lumi <= 252))
    return true;
  if (run == 208307 && (lumi >= 256 && lumi <= 259))
    return true;
  if (run == 208307 && (lumi >= 262 && lumi <= 275))
    return true;
  if (run == 208307 && (lumi >= 278 && lumi <= 342))
    return true;
  if (run == 208307 && (lumi >= 345 && lumi <= 450))
    return true;
  if (run == 208307 && (lumi >= 453 && lumi <= 527))
    return true;
  if (run == 208307 && (lumi >= 530 && lumi <= 583))
    return true;
  if (run == 208307 && (lumi >= 586 && lumi <= 605))
    return true;
  if (run == 208307 && (lumi >= 608 && lumi <= 616))
    return true;
  if (run == 208307 && (lumi >= 618 && lumi <= 667))
    return true;
  if (run == 208307 && (lumi >= 670 && lumi <= 761))
    return true;
  if (run == 208307 && (lumi >= 763 && lumi <= 798))
    return true;
  if (run == 208307 && (lumi >= 800 && lumi <= 889))
    return true;
  if (run == 208307 && (lumi >= 891 && lumi <= 893))
    return true;
  if (run == 208307 && (lumi >= 896 && lumi <= 1055))
    return true;
  if (run == 208307 && (lumi >= 1057 && lumi <= 1205))
    return true;
  if (run == 208307 && (lumi >= 1208 && lumi <= 1294))
    return true;
  if (run == 208307 && (lumi >= 1297 && lumi <= 1328))
    return true;
  if (run == 202013 && (lumi >= 1 && lumi <= 2))
    return true;
  if (run == 202013 && (lumi >= 5 && lumi <= 35))
    return true;
  if (run == 202013 && (lumi >= 38 && lumi <= 57))
    return true;
  if (run == 202012 && (lumi >= 98 && lumi <= 121))
    return true;
  if (run == 202012 && (lumi >= 126 && lumi <= 131))
    return true;
  if (run == 206187 && (lumi >= 107 && lumi <= 169))
    return true;
  if (run == 206187 && (lumi >= 172 && lumi <= 242))
    return true;
  if (run == 206187 && (lumi >= 245 && lumi <= 288))
    return true;
  if (run == 206187 && (lumi >= 290 && lumi <= 340))
    return true;
  if (run == 206187 && (lumi >= 343 && lumi <= 427))
    return true;
  if (run == 206187 && (lumi >= 429 && lumi <= 435))
    return true;
  if (run == 206187 && (lumi >= 437 && lumi <= 486))
    return true;
  if (run == 206187 && (lumi >= 489 && lumi <= 569))
    return true;
  if (run == 206187 && (lumi >= 571 && lumi <= 647))
    return true;
  if (run == 206187 && (lumi >= 649 && lumi <= 662))
    return true;
  if (run == 206187 && (lumi >= 664 && lumi <= 708))
    return true;
  if (run == 194153 && (lumi >= 1 && lumi <= 115))
    return true;
  if (run == 194150 && (lumi >= 42 && lumi <= 127))
    return true;
  if (run == 194150 && (lumi >= 129 && lumi <= 261))
    return true;
  if (run == 194150 && (lumi >= 264 && lumi <= 311))
    return true;
  if (run == 194151 && (lumi >= 47 && lumi <= 72))
    return true;
  if (run == 194151 && (lumi >= 75 && lumi <= 191))
    return true;
  if (run == 194151 && (lumi >= 193 && lumi <= 238))
    return true;
  if (run == 194151 && (lumi >= 240 && lumi <= 617))
    return true;
  if (run == 194151 && (lumi >= 619 && lumi <= 619))
    return true;
  if (run == 194151 && (lumi >= 621 && lumi <= 621))
    return true;
  if (run == 194151 && (lumi >= 623 && lumi <= 623))
    return true;
  if (run == 207214 && (lumi >= 57 && lumi <= 112))
    return true;
  if (run == 207214 && (lumi >= 114 && lumi <= 177))
    return true;
  if (run == 207214 && (lumi >= 179 && lumi <= 181))
    return true;
  if (run == 207214 && (lumi >= 184 && lumi <= 196))
    return true;
  if (run == 207214 && (lumi >= 199 && lumi <= 220))
    return true;
  if (run == 207214 && (lumi >= 223 && lumi <= 262))
    return true;
  if (run == 207214 && (lumi >= 265 && lumi <= 405))
    return true;
  if (run == 207214 && (lumi >= 408 && lumi <= 482))
    return true;
  if (run == 207214 && (lumi >= 485 && lumi <= 640))
    return true;
  if (run == 207214 && (lumi >= 643 && lumi <= 708))
    return true;
  if (run == 207214 && (lumi >= 718 && lumi <= 757))
    return true;
  if (run == 207214 && (lumi >= 759 && lumi <= 808))
    return true;
  if (run == 207214 && (lumi >= 811 && lumi <= 829))
    return true;
  if (run == 207217 && (lumi >= 1 && lumi <= 32))
    return true;
  if (run == 196453 && (lumi >= 1 && lumi <= 74))
    return true;
  if (run == 196453 && (lumi >= 77 && lumi <= 145))
    return true;
  if (run == 196453 && (lumi >= 147 && lumi <= 669))
    return true;
  if (run == 196453 && (lumi >= 673 && lumi <= 714))
    return true;
  if (run == 196453 && (lumi >= 717 && lumi <= 799))
    return true;
  if (run == 196453 && (lumi >= 802 && lumi <= 988))
    return true;
  if (run == 196453 && (lumi >= 991 && lumi <= 1178))
    return true;
  if (run == 196453 && (lumi >= 1180 && lumi <= 1180))
    return true;
  if (run == 196453 && (lumi >= 1182 && lumi <= 1248))
    return true;
  if (run == 196453 && (lumi >= 1250 && lumi <= 1528))
    return true;
  if (run == 196453 && (lumi >= 1531 && lumi <= 1647))
    return true;
  if (run == 207886 && (lumi >= 1 && lumi <= 30))
    return true;
  if (run == 207886 && (lumi >= 32 && lumi <= 90))
    return true;
  if (run == 207886 && (lumi >= 92 && lumi <= 156))
    return true;
  if (run == 207886 && (lumi >= 158 && lumi <= 166))
    return true;
  if (run == 207886 && (lumi >= 168 && lumi <= 171))
    return true;
  if (run == 207219 && (lumi >= 1 && lumi <= 112))
    return true;
  if (run == 195540 && (lumi >= 68 && lumi <= 123))
    return true;
  if (run == 195540 && (lumi >= 126 && lumi <= 137))
    return true;
  if (run == 195540 && (lumi >= 140 && lumi <= 283))
    return true;
  if (run == 195540 && (lumi >= 286 && lumi <= 319))
    return true;
  if (run == 206302 && (lumi >= 1 && lumi <= 8))
    return true;
  if (run == 206302 && (lumi >= 11 && lumi <= 33))
    return true;
  if (run == 206302 && (lumi >= 36 && lumi <= 44))
    return true;
  if (run == 206302 && (lumi >= 47 && lumi <= 82))
    return true;
  if (run == 206302 && (lumi >= 84 && lumi <= 108))
    return true;
  if (run == 206302 && (lumi >= 110 && lumi <= 149))
    return true;
  if (run == 206302 && (lumi >= 151 && lumi <= 186))
    return true;
  if (run == 206302 && (lumi >= 189 && lumi <= 229))
    return true;
  if (run == 206302 && (lumi >= 231 && lumi <= 232))
    return true;
  if (run == 206302 && (lumi >= 234 && lumi <= 241))
    return true;
  if (run == 206302 && (lumi >= 243 && lumi <= 276))
    return true;
  if (run == 195251 && (lumi >= 1 && lumi <= 131))
    return true;
  if (run == 195251 && (lumi >= 134 && lumi <= 137))
    return true;
  if (run == 195251 && (lumi >= 140 && lumi <= 152))
    return true;
  if (run == 195251 && (lumi >= 154 && lumi <= 165))
    return true;
  if (run == 195251 && (lumi >= 167 && lumi <= 242))
    return true;
  if (run == 205718 && (lumi >= 49 && lumi <= 75))
    return true;
  if (run == 205718 && (lumi >= 78 && lumi <= 97))
    return true;
  if (run == 205718 && (lumi >= 100 && lumi <= 103))
    return true;
  if (run == 205718 && (lumi >= 105 && lumi <= 176))
    return true;
  if (run == 205718 && (lumi >= 178 && lumi <= 338))
    return true;
  if (run == 205718 && (lumi >= 341 && lumi <= 361))
    return true;
  if (run == 205718 && (lumi >= 363 && lumi <= 524))
    return true;
  if (run == 205718 && (lumi >= 527 && lumi <= 531))
    return true;
  if (run == 205718 && (lumi >= 534 && lumi <= 589))
    return true;
  if (run == 205718 && (lumi >= 591 && lumi <= 694))
    return true;
  if (run == 205921 && (lumi >= 22 && lumi <= 73))
    return true;
  if (run == 205921 && (lumi >= 76 && lumi <= 268))
    return true;
  if (run == 205921 && (lumi >= 271 && lumi <= 394))
    return true;
  if (run == 205921 && (lumi >= 397 && lumi <= 401))
    return true;
  if (run == 205921 && (lumi >= 410 && lumi <= 428))
    return true;
  if (run == 205921 && (lumi >= 431 && lumi <= 498))
    return true;
  if (run == 205921 && (lumi >= 500 && lumi <= 571))
    return true;
  if (run == 205921 && (lumi >= 574 && lumi <= 779))
    return true;
  if (run == 205921 && (lumi >= 782 && lumi <= 853))
    return true;
  if (run == 190703 && (lumi >= 1 && lumi <= 252))
    return true;
  if (run == 190702 && (lumi >= 51 && lumi <= 53))
    return true;
  if (run == 190702 && (lumi >= 55 && lumi <= 122))
    return true;
  if (run == 190702 && (lumi >= 124 && lumi <= 169))
    return true;
  if (run == 190707 && (lumi >= 1 && lumi <= 237))
    return true;
  if (run == 190707 && (lumi >= 239 && lumi <= 257))
    return true;
  if (run == 190706 && (lumi >= 1 && lumi <= 126))
    return true;
  if (run == 199409 && (lumi >= 25 && lumi <= 54))
    return true;
  if (run == 199409 && (lumi >= 56 && lumi <= 89))
    return true;
  if (run == 199409 && (lumi >= 91 && lumi <= 204))
    return true;
  if (run == 199409 && (lumi >= 206 && lumi <= 290))
    return true;
  if (run == 199409 && (lumi >= 293 && lumi <= 583))
    return true;
  if (run == 199409 && (lumi >= 586 && lumi <= 602))
    return true;
  if (run == 199409 && (lumi >= 604 && lumi <= 1014))
    return true;
  if (run == 199409 && (lumi >= 1016 && lumi <= 1300))
    return true;
  if (run == 190704 && (lumi >= 1 && lumi <= 3))
    return true;
  if (run == 196364 && (lumi >= 1 && lumi <= 93))
    return true;
  if (run == 196364 && (lumi >= 96 && lumi <= 136))
    return true;
  if (run == 196364 && (lumi >= 139 && lumi <= 365))
    return true;
  if (run == 196364 && (lumi >= 368 && lumi <= 380))
    return true;
  if (run == 196364 && (lumi >= 382 && lumi <= 601))
    return true;
  if (run == 196364 && (lumi >= 603 && lumi <= 795))
    return true;
  if (run == 196364 && (lumi >= 798 && lumi <= 884))
    return true;
  if (run == 196364 && (lumi >= 887 && lumi <= 1196))
    return true;
  if (run == 196364 && (lumi >= 1199 && lumi <= 1200))
    return true;
  if (run == 196364 && (lumi >= 1203 && lumi <= 1299))
    return true;
  if (run == 206866 && (lumi >= 1 && lumi <= 30))
    return true;
  if (run == 206866 && (lumi >= 33 && lumi <= 113))
    return true;
  if (run == 206866 && (lumi >= 115 && lumi <= 274))
    return true;
  if (run == 194303 && (lumi >= 56 && lumi <= 66))
    return true;
  if (run == 194303 && (lumi >= 69 && lumi <= 102))
    return true;
  if (run == 196218 && (lumi >= 55 && lumi <= 199))
    return true;
  if (run == 196218 && (lumi >= 201 && lumi <= 224))
    return true;
  if (run == 196218 && (lumi >= 226 && lumi <= 393))
    return true;
  if (run == 196218 && (lumi >= 396 && lumi <= 494))
    return true;
  if (run == 196218 && (lumi >= 496 && lumi <= 741))
    return true;
  if (run == 196218 && (lumi >= 744 && lumi <= 752))
    return true;
  if (run == 196218 && (lumi >= 754 && lumi <= 757))
    return true;
  if (run == 196218 && (lumi >= 759 && lumi <= 820))
    return true;
  if (run == 194305 && (lumi >= 1 && lumi <= 84))
    return true;
  if (run == 194304 && (lumi >= 1 && lumi <= 43))
    return true;
  if (run == 194304 && (lumi >= 46 && lumi <= 46))
    return true;
  if (run == 193336 && (lumi >= 1 && lumi <= 264))
    return true;
  if (run == 193336 && (lumi >= 267 && lumi <= 492))
    return true;
  if (run == 193336 && (lumi >= 495 && lumi <= 684))
    return true;
  if (run == 193336 && (lumi >= 687 && lumi <= 729))
    return true;
  if (run == 193336 && (lumi >= 732 && lumi <= 951))
    return true;
  if (run == 201669 && (lumi >= 1 && lumi <= 9))
    return true;
  if (run == 201669 && (lumi >= 12 && lumi <= 136))
    return true;
  if (run == 201669 && (lumi >= 139 && lumi <= 141))
    return true;
  if (run == 201669 && (lumi >= 143 && lumi <= 165))
    return true;
  if (run == 201668 && (lumi >= 78 && lumi <= 157))
    return true;
  if (run == 195099 && (lumi >= 70 && lumi <= 144))
    return true;
  if (run == 195099 && (lumi >= 147 && lumi <= 186))
    return true;
  if (run == 195099 && (lumi >= 189 && lumi <= 208))
    return true;
  if (run == 195099 && (lumi >= 211 && lumi <= 224))
    return true;
  if (run == 195099 && (lumi >= 227 && lumi <= 248))
    return true;
  if (run == 206466 && (lumi >= 24 && lumi <= 137))
    return true;
  if (run == 206466 && (lumi >= 140 && lumi <= 277))
    return true;
  if (run == 206466 && (lumi >= 280 && lumi <= 296))
    return true;
  if (run == 206466 && (lumi >= 299 && lumi <= 303))
    return true;
  if (run == 206466 && (lumi >= 306 && lumi <= 405))
    return true;
  if (run == 206466 && (lumi >= 407 && lumi <= 419))
    return true;
  if (run == 206466 && (lumi >= 422 && lumi <= 477))
    return true;
  if (run == 206466 && (lumi >= 480 && lumi <= 511))
    return true;
  if (run == 206466 && (lumi >= 514 && lumi <= 676))
    return true;
  if (run == 190906 && (lumi >= 73 && lumi <= 256))
    return true;
  if (run == 190906 && (lumi >= 259 && lumi <= 354))
    return true;
  if (run == 190906 && (lumi >= 356 && lumi <= 496))
    return true;
  if (run == 200041 && (lumi >= 62 && lumi <= 83))
    return true;
  if (run == 200041 && (lumi >= 85 && lumi <= 157))
    return true;
  if (run == 200041 && (lumi >= 162 && lumi <= 274))
    return true;
  if (run == 200041 && (lumi >= 277 && lumi <= 318))
    return true;
  if (run == 200041 && (lumi >= 321 && lumi <= 335))
    return true;
  if (run == 200041 && (lumi >= 337 && lumi <= 386))
    return true;
  if (run == 200041 && (lumi >= 388 && lumi <= 389))
    return true;
  if (run == 200041 && (lumi >= 392 && lumi <= 400))
    return true;
  if (run == 200041 && (lumi >= 402 && lumi <= 568))
    return true;
  if (run == 200041 && (lumi >= 571 && lumi <= 593))
    return true;
  if (run == 200041 && (lumi >= 595 && lumi <= 646))
    return true;
  if (run == 200041 && (lumi >= 649 && lumi <= 728))
    return true;
  if (run == 200041 && (lumi >= 731 && lumi <= 860))
    return true;
  if (run == 200041 && (lumi >= 862 && lumi <= 930))
    return true;
  if (run == 200041 && (lumi >= 932 && lumi <= 1096))
    return true;
  if (run == 200042 && (lumi >= 1 && lumi <= 110))
    return true;
  if (run == 200042 && (lumi >= 112 && lumi <= 536))
    return true;
  if (run == 204567 && (lumi >= 1 && lumi <= 38))
    return true;
  if (run == 204566 && (lumi >= 1 && lumi <= 12))
    return true;
  if (run == 204565 && (lumi >= 1 && lumi <= 48))
    return true;
  if (run == 203709 && (lumi >= 1 && lumi <= 121))
    return true;
  if (run == 204563 && (lumi >= 91 && lumi <= 99))
    return true;
  if (run == 204563 && (lumi >= 102 && lumi <= 178))
    return true;
  if (run == 204563 && (lumi >= 180 && lumi <= 219))
    return true;
  if (run == 204563 && (lumi >= 222 && lumi <= 229))
    return true;
  if (run == 204563 && (lumi >= 231 && lumi <= 364))
    return true;
  if (run == 204563 && (lumi >= 366 && lumi <= 366))
    return true;
  if (run == 204563 && (lumi >= 369 && lumi <= 470))
    return true;
  if (run == 204563 && (lumi >= 473 && lumi <= 524))
    return true;
  if (run == 204563 && (lumi >= 527 && lumi <= 571))
    return true;
  if (run == 194051 && (lumi >= 1 && lumi <= 12))
    return true;
  if (run == 206484 && (lumi >= 73 && lumi <= 95))
    return true;
  if (run == 206484 && (lumi >= 98 && lumi <= 133))
    return true;
  if (run == 206484 && (lumi >= 136 && lumi <= 163))
    return true;
  if (run == 206484 && (lumi >= 166 && lumi <= 186))
    return true;
  if (run == 206484 && (lumi >= 189 && lumi <= 384))
    return true;
  if (run == 206484 && (lumi >= 387 && lumi <= 463))
    return true;
  if (run == 206484 && (lumi >= 465 && lumi <= 551))
    return true;
  if (run == 206484 && (lumi >= 554 && lumi <= 554))
    return true;
  if (run == 206484 && (lumi >= 556 && lumi <= 669))
    return true;
  if (run == 194119 && (lumi >= 1 && lumi <= 229))
    return true;
  if (run == 194119 && (lumi >= 232 && lumi <= 261))
    return true;
  if (run == 198049 && (lumi >= 1 && lumi <= 11))
    return true;
  if (run == 198049 && (lumi >= 14 && lumi <= 57))
    return true;
  if (run == 202333 && (lumi >= 1 && lumi <= 235))
    return true;
  if (run == 195757 && (lumi >= 1 && lumi <= 82))
    return true;
  if (run == 195757 && (lumi >= 85 && lumi <= 115))
    return true;
  if (run == 195757 && (lumi >= 118 && lumi <= 161))
    return true;
  if (run == 195757 && (lumi >= 163 && lumi <= 206))
    return true;
  if (run == 190679 && (lumi >= 1 && lumi <= 55))
    return true;
  if (run == 195758 && (lumi >= 1 && lumi <= 18))
    return true;
  if (run == 195164 && (lumi >= 1 && lumi <= 64))
    return true;
  if (run == 195165 && (lumi >= 1 && lumi <= 4))
    return true;
  if (run == 195165 && (lumi >= 7 && lumi <= 41))
    return true;
  if (run == 195165 && (lumi >= 44 && lumi <= 54))
    return true;
  if (run == 195165 && (lumi >= 56 && lumi <= 153))
    return true;
  if (run == 195165 && (lumi >= 156 && lumi <= 260))
    return true;
  if (run == 195165 && (lumi >= 263 && lumi <= 266))
    return true;
  if (run == 201062 && (lumi >= 78 && lumi <= 268))
    return true;
  if (run == 195163 && (lumi >= 72 && lumi <= 138))
    return true;
  if (run == 195163 && (lumi >= 140 && lumi <= 224))
    return true;
  if (run == 195163 && (lumi >= 227 && lumi <= 240))
    return true;
  if (run == 195163 && (lumi >= 243 && lumi <= 243))
    return true;
  if (run == 195163 && (lumi >= 246 && lumi <= 347))
    return true;
  if (run == 194270 && (lumi >= 56 && lumi <= 68))
    return true;
  if (run == 199564 && (lumi >= 1 && lumi <= 3))
    return true;
  if (run == 203832 && (lumi >= 1 && lumi <= 11))
    return true;
  if (run == 203833 && (lumi >= 1 && lumi <= 70))
    return true;
  if (run == 203833 && (lumi >= 73 && lumi <= 128))
    return true;
  if (run == 203834 && (lumi >= 1 && lumi <= 40))
    return true;
  if (run == 203835 && (lumi >= 1 && lumi <= 70))
    return true;
  if (run == 203835 && (lumi >= 73 && lumi <= 358))
    return true;
  if (run == 201692 && (lumi >= 78 && lumi <= 81))
    return true;
  if (run == 201692 && (lumi >= 83 && lumi <= 179))
    return true;
  if (run == 199569 && (lumi >= 1 && lumi <= 2))
    return true;
  if (run == 199569 && (lumi >= 5 && lumi <= 136))
    return true;
  if (run == 199569 && (lumi >= 139 && lumi <= 367))
    return true;
  if (run == 206304 && (lumi >= 1 && lumi <= 4))
    return true;
  if (run == 206304 && (lumi >= 6 && lumi <= 62))
    return true;
  if (run == 198215 && (lumi >= 1 && lumi <= 12))
    return true;
  if (run == 202205 && (lumi >= 94 && lumi <= 114))
    return true;
  if (run == 198210 && (lumi >= 1 && lumi <= 221))
    return true;
  if (run == 198213 && (lumi >= 1 && lumi <= 107))
    return true;
  if (run == 198212 && (lumi >= 1 && lumi <= 574))
    return true;
  if (run == 191271 && (lumi >= 56 && lumi <= 223))
    return true;
  if (run == 191271 && (lumi >= 225 && lumi <= 363))
    return true;
  if (run == 205193 && (lumi >= 80 && lumi <= 109))
    return true;
  if (run == 205193 && (lumi >= 111 && lumi <= 349))
    return true;
  if (run == 205193 && (lumi >= 352 && lumi <= 486))
    return true;
  if (run == 205193 && (lumi >= 488 && lumi <= 650))
    return true;
  if (run == 205193 && (lumi >= 652 && lumi <= 712))
    return true;
  if (run == 205193 && (lumi >= 714 && lumi <= 902))
    return true;
  if (run == 191720 && (lumi >= 1 && lumi <= 1))
    return true;
  if (run == 191720 && (lumi >= 3 && lumi <= 15))
    return true;
  if (run == 191720 && (lumi >= 17 && lumi <= 181))
    return true;
  if (run == 191276 && (lumi >= 1 && lumi <= 16))
    return true;
  if (run == 191277 && (lumi >= 1 && lumi <= 28))
    return true;
  if (run == 191277 && (lumi >= 30 && lumi <= 164))
    return true;
  if (run == 191277 && (lumi >= 167 && lumi <= 253))
    return true;
  if (run == 191277 && (lumi >= 255 && lumi <= 457))
    return true;
  if (run == 191277 && (lumi >= 460 && lumi <= 535))
    return true;
  if (run == 191277 && (lumi >= 537 && lumi <= 576))
    return true;
  if (run == 191277 && (lumi >= 579 && lumi <= 775))
    return true;
  if (run == 191277 && (lumi >= 778 && lumi <= 811))
    return true;
  if (run == 191277 && (lumi >= 813 && lumi <= 849))
    return true;
  if (run == 191833 && (lumi >= 1 && lumi <= 1))
    return true;
  if (run == 191833 && (lumi >= 3 && lumi <= 103))
    return true;
  if (run == 191830 && (lumi >= 54 && lumi <= 242))
    return true;
  if (run == 191830 && (lumi >= 245 && lumi <= 301))
    return true;
  if (run == 191830 && (lumi >= 304 && lumi <= 393))
    return true;
  if (run == 191837 && (lumi >= 1 && lumi <= 44))
    return true;
  if (run == 191837 && (lumi >= 47 && lumi <= 53))
    return true;
  if (run == 191837 && (lumi >= 56 && lumi <= 65))
    return true;
  if (run == 191834 && (lumi >= 1 && lumi <= 30))
    return true;
  if (run == 191834 && (lumi >= 33 && lumi <= 74))
    return true;
  if (run == 191834 && (lumi >= 77 && lumi <= 299))
    return true;
  if (run == 191834 && (lumi >= 302 && lumi <= 352))
    return true;
  if (run == 207905 && (lumi >= 75 && lumi <= 196))
    return true;
  if (run == 207905 && (lumi >= 198 && lumi <= 281))
    return true;
  if (run == 207905 && (lumi >= 284 && lumi <= 329))
    return true;
  if (run == 207905 && (lumi >= 331 && lumi <= 402))
    return true;
  if (run == 207905 && (lumi >= 404 && lumi <= 565))
    return true;
  if (run == 207905 && (lumi >= 568 && lumi <= 672))
    return true;
  if (run == 207905 && (lumi >= 675 && lumi <= 805))
    return true;
  if (run == 207905 && (lumi >= 807 && lumi <= 850))
    return true;
  if (run == 207905 && (lumi >= 852 && lumi <= 861))
    return true;
  if (run == 207905 && (lumi >= 864 && lumi <= 884))
    return true;
  if (run == 207905 && (lumi >= 886 && lumi <= 1180))
    return true;
  if (run == 207905 && (lumi >= 1183 && lumi <= 1283))
    return true;
  if (run == 207905 && (lumi >= 1285 && lumi <= 1331))
    return true;
  if (run == 207905 && (lumi >= 1333 && lumi <= 1515))
    return true;
  if (run == 207905 && (lumi >= 1518 && lumi <= 1734))
    return true;
  if (run == 207905 && (lumi >= 1737 && lumi <= 1796))
    return true;
  if (run == 206188 && (lumi >= 1 && lumi <= 40))
    return true;
  if (run == 206188 && (lumi >= 42 && lumi <= 55))
    return true;
  if (run == 205667 && (lumi >= 1 && lumi <= 165))
    return true;
  if (run == 205667 && (lumi >= 168 && lumi <= 282))
    return true;
  if (run == 205667 && (lumi >= 285 && lumi <= 318))
    return true;
  if (run == 205667 && (lumi >= 321 && lumi <= 412))
    return true;
  if (run == 205667 && (lumi >= 415 && lumi <= 689))
    return true;
  if (run == 205667 && (lumi >= 692 && lumi <= 751))
    return true;
  if (run == 205667 && (lumi >= 754 && lumi <= 774))
    return true;
  if (run == 205667 && (lumi >= 777 && lumi <= 1109))
    return true;
  if (run == 205666 && (lumi >= 60 && lumi <= 119))
    return true;
  if (run == 205666 && (lumi >= 122 && lumi <= 165))
    return true;
  if (run == 205666 && (lumi >= 168 && lumi <= 259))
    return true;
  if (run == 205666 && (lumi >= 261 && lumi <= 322))
    return true;
  if (run == 205666 && (lumi >= 325 && lumi <= 578))
    return true;
  if (run == 205666 && (lumi >= 580 && lumi <= 594))
    return true;
  if (run == 205666 && (lumi >= 597 && lumi <= 721))
    return true;
  if (run == 205666 && (lumi >= 724 && lumi <= 739))
    return true;
  if (run == 201115 && (lumi >= 1 && lumi <= 73))
    return true;
  if (run == 201114 && (lumi >= 1 && lumi <= 14))
    return true;
  if (run == 194117 && (lumi >= 1 && lumi <= 38))
    return true;
  if (run == 205683 && (lumi >= 76 && lumi <= 82))
    return true;
  if (run == 205683 && (lumi >= 85 && lumi <= 178))
    return true;
  if (run == 205683 && (lumi >= 181 && lumi <= 198))
    return true;
  if (run == 205683 && (lumi >= 201 && lumi <= 305))
    return true;
  if (run == 194115 && (lumi >= 66 && lumi <= 184))
    return true;
  if (run == 194115 && (lumi >= 186 && lumi <= 338))
    return true;
  if (run == 194115 && (lumi >= 340 && lumi <= 346))
    return true;
  if (run == 194115 && (lumi >= 348 && lumi <= 493))
    return true;
  if (run == 194115 && (lumi >= 496 && lumi <= 731))
    return true;
  if (run == 194115 && (lumi >= 819 && lumi <= 857))
    return true;
  if (run == 194052 && (lumi >= 1 && lumi <= 99))
    return true;
  if (run == 194052 && (lumi >= 102 && lumi <= 166))
    return true;
  if (run == 200532 && (lumi >= 1 && lumi <= 37))
    return true;
  if (run == 194050 && (lumi >= 53 && lumi <= 113))
    return true;
  if (run == 194050 && (lumi >= 116 && lumi <= 273))
    return true;
  if (run == 194050 && (lumi >= 275 && lumi <= 355))
    return true;
  if (run == 194050 && (lumi >= 357 && lumi <= 369))
    return true;
  if (run == 194050 && (lumi >= 372 && lumi <= 391))
    return true;
  if (run == 194050 && (lumi >= 394 && lumi <= 490))
    return true;
  if (run == 194050 && (lumi >= 492 && lumi <= 814))
    return true;
  if (run == 194050 && (lumi >= 816 && lumi <= 1435))
    return true;
  if (run == 194050 && (lumi >= 1437 && lumi <= 1735))
    return true;
  if (run == 194050 && (lumi >= 1760 && lumi <= 1888))
    return true;
  if (run == 206303 && (lumi >= 1 && lumi <= 19))
    return true;
  if (run == 206303 && (lumi >= 23 && lumi <= 286))
    return true;
  if (run == 194704 && (lumi >= 1 && lumi <= 41))
    return true;
  if (run == 194704 && (lumi >= 44 && lumi <= 545))
    return true;
  if (run == 194704 && (lumi >= 548 && lumi <= 592))
    return true;
  if (run == 196437 && (lumi >= 1 && lumi <= 1))
    return true;
  if (run == 196437 && (lumi >= 3 && lumi <= 74))
    return true;
  if (run == 196437 && (lumi >= 77 && lumi <= 169))
    return true;
  if (run == 194702 && (lumi >= 1 && lumi <= 138))
    return true;
  if (run == 194702 && (lumi >= 141 && lumi <= 191))
    return true;
  if (run == 202116 && (lumi >= 59 && lumi <= 60))
    return true;
  if (run == 207398 && (lumi >= 1 && lumi <= 14))
    return true;
  if (run == 207398 && (lumi >= 16 && lumi <= 33))
    return true;
  if (run == 196438 && (lumi >= 1 && lumi <= 181))
    return true;
  if (run == 196438 && (lumi >= 184 && lumi <= 699))
    return true;
  if (run == 196438 && (lumi >= 701 && lumi <= 1269))
    return true;
  if (run == 202054 && (lumi >= 6 && lumi <= 266))
    return true;
  if (run == 202054 && (lumi >= 268 && lumi <= 489))
    return true;
  if (run == 202054 && (lumi >= 492 && lumi <= 605))
    return true;
  if (run == 202054 && (lumi >= 608 && lumi <= 631))
    return true;
  if (run == 203909 && (lumi >= 79 && lumi <= 113))
    return true;
  if (run == 203909 && (lumi >= 116 && lumi <= 117))
    return true;
  if (run == 203909 && (lumi >= 120 && lumi <= 140))
    return true;
  if (run == 203909 && (lumi >= 143 && lumi <= 382))
    return true;
  if (run == 194464 && (lumi >= 1 && lumi <= 127))
    return true;
  if (run == 194464 && (lumi >= 130 && lumi <= 142))
    return true;
  if (run == 194464 && (lumi >= 145 && lumi <= 210))
    return true;
  if (run == 196531 && (lumi >= 62 && lumi <= 150))
    return true;
  if (run == 196531 && (lumi >= 152 && lumi <= 253))
    return true;
  if (run == 196531 && (lumi >= 256 && lumi <= 285))
    return true;
  if (run == 196531 && (lumi >= 288 && lumi <= 302))
    return true;
  if (run == 196531 && (lumi >= 305 && lumi <= 422))
    return true;
  if (run == 196531 && (lumi >= 425 && lumi <= 440))
    return true;
  if (run == 207099 && (lumi >= 83 && lumi <= 134))
    return true;
  if (run == 207099 && (lumi >= 137 && lumi <= 172))
    return true;
  if (run == 207099 && (lumi >= 175 && lumi <= 213))
    return true;
  if (run == 207099 && (lumi >= 216 && lumi <= 314))
    return true;
  if (run == 207099 && (lumi >= 316 && lumi <= 320))
    return true;
  if (run == 207099 && (lumi >= 323 && lumi <= 330))
    return true;
  if (run == 207099 && (lumi >= 333 && lumi <= 367))
    return true;
  if (run == 207099 && (lumi >= 370 && lumi <= 481))
    return true;
  if (run == 207099 && (lumi >= 484 && lumi <= 602))
    return true;
  if (run == 207099 && (lumi >= 605 && lumi <= 755))
    return true;
  if (run == 207099 && (lumi >= 757 && lumi <= 1046))
    return true;
  if (run == 207099 && (lumi >= 1048 && lumi <= 1171))
    return true;
  if (run == 200229 && (lumi >= 1 && lumi <= 33))
    return true;
  if (run == 200229 && (lumi >= 41 && lumi <= 219))
    return true;
  if (run == 200229 && (lumi >= 222 && lumi <= 244))
    return true;
  if (run == 200229 && (lumi >= 247 && lumi <= 290))
    return true;
  if (run == 200229 && (lumi >= 293 && lumi <= 624))
    return true;
  if (run == 200229 && (lumi >= 627 && lumi <= 629))
    return true;
  if (run == 205618 && (lumi >= 1 && lumi <= 12))
    return true;
  if (run == 205617 && (lumi >= 1 && lumi <= 29))
    return true;
  if (run == 205617 && (lumi >= 32 && lumi <= 102))
    return true;
  if (run == 205617 && (lumi >= 105 && lumi <= 123))
    return true;
  if (run == 205617 && (lumi >= 125 && lumi <= 140))
    return true;
  if (run == 205617 && (lumi >= 143 && lumi <= 264))
    return true;
  if (run == 205617 && (lumi >= 266 && lumi <= 448))
    return true;
  if (run == 205617 && (lumi >= 451 && lumi <= 532))
    return true;
  if (run == 205617 && (lumi >= 534 && lumi <= 547))
    return true;
  if (run == 205614 && (lumi >= 1 && lumi <= 4))
    return true;
  if (run == 205614 && (lumi >= 7 && lumi <= 40))
    return true;
  if (run == 200152 && (lumi >= 74 && lumi <= 116))
    return true;
  if (run == 195947 && (lumi >= 23 && lumi <= 62))
    return true;
  if (run == 195947 && (lumi >= 64 && lumi <= 88))
    return true;
  if (run == 195948 && (lumi >= 51 && lumi <= 116))
    return true;
  if (run == 195948 && (lumi >= 119 && lumi <= 144))
    return true;
  if (run == 195948 && (lumi >= 147 && lumi <= 147))
    return true;
  if (run == 195948 && (lumi >= 150 && lumi <= 352))
    return true;
  if (run == 195948 && (lumi >= 355 && lumi <= 369))
    return true;
  if (run == 195948 && (lumi >= 372 && lumi <= 402))
    return true;
  if (run == 195948 && (lumi >= 404 && lumi <= 500))
    return true;
  if (run == 195948 && (lumi >= 503 && lumi <= 540))
    return true;
  if (run == 195948 && (lumi >= 543 && lumi <= 565))
    return true;
  if (run == 195948 && (lumi >= 567 && lumi <= 602))
    return true;
  if (run == 195948 && (lumi >= 605 && lumi <= 615))
    return true;
  if (run == 191695 && (lumi >= 1 && lumi <= 1))
    return true;
  if (run == 191718 && (lumi >= 43 && lumi <= 95))
    return true;
  if (run == 191718 && (lumi >= 98 && lumi <= 207))
    return true;
  if (run == 190945 && (lumi >= 124 && lumi <= 207))
    return true;
  if (run == 191056 && (lumi >= 1 && lumi <= 1))
    return true;
  if (run == 191056 && (lumi >= 4 && lumi <= 9))
    return true;
  if (run == 191056 && (lumi >= 16 && lumi <= 17))
    return true;
  if (run == 191056 && (lumi >= 19 && lumi <= 19))
    return true;
  if (run == 191057 && (lumi >= 1 && lumi <= 1))
    return true;
  if (run == 191057 && (lumi >= 4 && lumi <= 40))
    return true;
  if (run == 190949 && (lumi >= 1 && lumi <= 81))
    return true;
  if (run == 201191 && (lumi >= 75 && lumi <= 98))
    return true;
  if (run == 201191 && (lumi >= 100 && lumi <= 216))
    return true;
  if (run == 201191 && (lumi >= 218 && lumi <= 389))
    return true;
  if (run == 201191 && (lumi >= 392 && lumi <= 492))
    return true;
  if (run == 201191 && (lumi >= 494 && lumi <= 506))
    return true;
  if (run == 201191 && (lumi >= 509 && lumi <= 585))
    return true;
  if (run == 201191 && (lumi >= 587 && lumi <= 594))
    return true;
  if (run == 201191 && (lumi >= 597 && lumi <= 607))
    return true;
  if (run == 201191 && (lumi >= 609 && lumi <= 794))
    return true;
  if (run == 201191 && (lumi >= 796 && lumi <= 838))
    return true;
  if (run == 201191 && (lumi >= 841 && lumi <= 974))
    return true;
  if (run == 201191 && (lumi >= 977 && lumi <= 1105))
    return true;
  if (run == 201191 && (lumi >= 1108 && lumi <= 1117))
    return true;
  if (run == 201191 && (lumi >= 1120 && lumi <= 1382))
    return true;
  if (run == 201191 && (lumi >= 1385 && lumi <= 1386))
    return true;
  if (run == 201625 && (lumi >= 211 && lumi <= 312))
    return true;
  if (run == 201625 && (lumi >= 315 && lumi <= 348))
    return true;
  if (run == 201625 && (lumi >= 351 && lumi <= 416))
    return true;
  if (run == 201625 && (lumi >= 418 && lumi <= 588))
    return true;
  if (run == 201625 && (lumi >= 591 && lumi <= 671))
    return true;
  if (run == 201625 && (lumi >= 673 && lumi <= 758))
    return true;
  if (run == 201625 && (lumi >= 760 && lumi <= 791))
    return true;
  if (run == 201625 && (lumi >= 793 && lumi <= 944))
    return true;
  if (run == 201624 && (lumi >= 83 && lumi <= 92))
    return true;
  if (run == 201624 && (lumi >= 95 && lumi <= 240))
    return true;
  if (run == 201624 && (lumi >= 270 && lumi <= 270))
    return true;
  if (run == 193124 && (lumi >= 1 && lumi <= 52))
    return true;
  if (run == 208339 && (lumi >= 77 && lumi <= 89))
    return true;
  if (run == 208339 && (lumi >= 91 && lumi <= 122))
    return true;
  if (run == 208339 && (lumi >= 125 && lumi <= 208))
    return true;
  if (run == 208339 && (lumi >= 211 && lumi <= 346))
    return true;
  if (run == 208339 && (lumi >= 349 && lumi <= 363))
    return true;
  if (run == 207320 && (lumi >= 1 && lumi <= 110))
    return true;
  if (run == 207320 && (lumi >= 112 && lumi <= 350))
    return true;
  if (run == 193123 && (lumi >= 1 && lumi <= 27))
    return true;
  if (run == 207898 && (lumi >= 1 && lumi <= 33))
    return true;
  if (run == 207898 && (lumi >= 36 && lumi <= 57))
    return true;
  if (run == 207898 && (lumi >= 60 && lumi <= 235))
    return true;
  if (run == 207898 && (lumi >= 239 && lumi <= 257))
    return true;
  if (run == 207898 && (lumi >= 260 && lumi <= 277))
    return true;
  if (run == 206542 && (lumi >= 1 && lumi <= 115))
    return true;
  if (run == 206542 && (lumi >= 117 && lumi <= 165))
    return true;
  if (run == 206542 && (lumi >= 168 && lumi <= 511))
    return true;
  if (run == 206542 && (lumi >= 514 && lumi <= 547))
    return true;
  if (run == 206542 && (lumi >= 550 && lumi <= 603))
    return true;
  if (run == 206542 && (lumi >= 606 && lumi <= 668))
    return true;
  if (run == 206542 && (lumi >= 671 && lumi <= 727))
    return true;
  if (run == 206542 && (lumi >= 730 && lumi <= 739))
    return true;
  if (run == 206542 && (lumi >= 741 && lumi <= 833))
    return true;
  if (run == 207487 && (lumi >= 50 && lumi <= 98))
    return true;
  if (run == 207487 && (lumi >= 101 && lumi <= 311))
    return true;
  if (run == 207487 && (lumi >= 313 && lumi <= 359))
    return true;
  if (run == 207487 && (lumi >= 363 && lumi <= 468))
    return true;
  if (run == 207487 && (lumi >= 471 && lumi <= 472))
    return true;
  if (run == 208686 && (lumi >= 73 && lumi <= 79))
    return true;
  if (run == 208686 && (lumi >= 82 && lumi <= 181))
    return true;
  if (run == 208686 && (lumi >= 183 && lumi <= 224))
    return true;
  if (run == 208686 && (lumi >= 227 && lumi <= 243))
    return true;
  if (run == 208686 && (lumi >= 246 && lumi <= 311))
    return true;
  if (run == 208686 && (lumi >= 313 && lumi <= 459))
    return true;
  if (run == 203992 && (lumi >= 1 && lumi <= 15))
    return true;
  if (run == 203994 && (lumi >= 1 && lumi <= 56))
    return true;
  if (run == 203994 && (lumi >= 59 && lumi <= 136))
    return true;
  if (run == 203994 && (lumi >= 139 && lumi <= 304))
    return true;
  if (run == 203994 && (lumi >= 306 && lumi <= 342))
    return true;
  if (run == 203994 && (lumi >= 344 && lumi <= 425))
    return true;
  if (run == 205339 && (lumi >= 71 && lumi <= 175))
    return true;
  if (run == 205339 && (lumi >= 178 && lumi <= 213))
    return true;
  if (run == 205339 && (lumi >= 216 && lumi <= 230))
    return true;
  if (run == 205339 && (lumi >= 233 && lumi <= 262))
    return true;
  if (run == 205339 && (lumi >= 265 && lumi <= 404))
    return true;
  if (run == 199436 && (lumi >= 1 && lumi <= 113))
    return true;
  if (run == 199436 && (lumi >= 116 && lumi <= 254))
    return true;
  if (run == 199436 && (lumi >= 257 && lumi <= 675))
    return true;
  if (run == 199436 && (lumi >= 678 && lumi <= 748))
    return true;
  if (run == 199435 && (lumi >= 63 && lumi <= 106))
    return true;
  if (run == 199435 && (lumi >= 109 && lumi <= 261))
    return true;
  if (run == 199435 && (lumi >= 263 && lumi <= 579))
    return true;
  if (run == 199435 && (lumi >= 582 && lumi <= 654))
    return true;
  if (run == 199435 && (lumi >= 656 && lumi <= 696))
    return true;
  if (run == 199435 && (lumi >= 699 && lumi <= 1034))
    return true;
  if (run == 199435 && (lumi >= 1037 && lumi <= 1144))
    return true;
  if (run == 199435 && (lumi >= 1147 && lumi <= 1327))
    return true;
  if (run == 199435 && (lumi >= 1330 && lumi <= 1411))
    return true;
  if (run == 199435 && (lumi >= 1414 && lumi <= 1431))
    return true;
  if (run == 199435 && (lumi >= 1434 && lumi <= 1441))
    return true;
  if (run == 199435 && (lumi >= 1444 && lumi <= 1487))
    return true;
  if (run == 199435 && (lumi >= 1489 && lumi <= 1610))
    return true;
  if (run == 204250 && (lumi >= 92 && lumi <= 118))
    return true;
  if (run == 204250 && (lumi >= 121 && lumi <= 177))
    return true;
  if (run == 204250 && (lumi >= 179 && lumi <= 285))
    return true;
  if (run == 204250 && (lumi >= 287 && lumi <= 336))
    return true;
  if (run == 204250 && (lumi >= 339 && lumi <= 400))
    return true;
  if (run == 204250 && (lumi >= 403 && lumi <= 521))
    return true;
  if (run == 204250 && (lumi >= 524 && lumi <= 543))
    return true;
  if (run == 204250 && (lumi >= 546 && lumi <= 682))
    return true;
  if (run == 204250 && (lumi >= 684 && lumi <= 801))
    return true;
  if (run == 190738 && (lumi >= 1 && lumi <= 130))
    return true;
  if (run == 190738 && (lumi >= 133 && lumi <= 226))
    return true;
  if (run == 190738 && (lumi >= 229 && lumi <= 349))
    return true;
  if (run == 190733 && (lumi >= 71 && lumi <= 96))
    return true;
  if (run == 190733 && (lumi >= 99 && lumi <= 389))
    return true;
  if (run == 190733 && (lumi >= 392 && lumi <= 460))
    return true;
  if (run == 190736 && (lumi >= 1 && lumi <= 80))
    return true;
  if (run == 190736 && (lumi >= 83 && lumi <= 185))
    return true;
  if (run == 201718 && (lumi >= 58 && lumi <= 108))
    return true;
  if (run == 199960 && (lumi >= 72 && lumi <= 139))
    return true;
  if (run == 199960 && (lumi >= 141 && lumi <= 197))
    return true;
  if (run == 199960 && (lumi >= 204 && lumi <= 232))
    return true;
  if (run == 199960 && (lumi >= 235 && lumi <= 363))
    return true;
  if (run == 199960 && (lumi >= 365 && lumi <= 367))
    return true;
  if (run == 199960 && (lumi >= 370 && lumi <= 380))
    return true;
  if (run == 199960 && (lumi >= 383 && lumi <= 459))
    return true;
  if (run == 199960 && (lumi >= 461 && lumi <= 466))
    return true;
  if (run == 199960 && (lumi >= 469 && lumi <= 485))
    return true;
  if (run == 198522 && (lumi >= 65 && lumi <= 144))
    return true;
  if (run == 198522 && (lumi >= 147 && lumi <= 208))
    return true;
  if (run == 201658 && (lumi >= 1 && lumi <= 19))
    return true;
  if (run == 201658 && (lumi >= 21 && lumi <= 118))
    return true;
  if (run == 201658 && (lumi >= 121 && lumi <= 136))
    return true;
  if (run == 201658 && (lumi >= 139 && lumi <= 288))
    return true;
  if (run == 199967 && (lumi >= 60 && lumi <= 120))
    return true;
  if (run == 199967 && (lumi >= 122 && lumi <= 170))
    return true;
  if (run == 199967 && (lumi >= 172 && lumi <= 198))
    return true;
  if (run == 201657 && (lumi >= 77 && lumi <= 93))
    return true;
  if (run == 201657 && (lumi >= 95 && lumi <= 108))
    return true;
  if (run == 201657 && (lumi >= 110 && lumi <= 118))
    return true;
  if (run == 198485 && (lumi >= 68 && lumi <= 109))
    return true;
  if (run == 198485 && (lumi >= 112 && lumi <= 134))
    return true;
  if (run == 198485 && (lumi >= 136 && lumi <= 181))
    return true;
  if (run == 198485 && (lumi >= 184 && lumi <= 239))
    return true;
  if (run == 198487 && (lumi >= 1 && lumi <= 145))
    return true;
  if (run == 198487 && (lumi >= 147 && lumi <= 514))
    return true;
  if (run == 198487 && (lumi >= 517 && lumi <= 668))
    return true;
  if (run == 198487 && (lumi >= 671 && lumi <= 733))
    return true;
  if (run == 198487 && (lumi >= 736 && lumi <= 757))
    return true;
  if (run == 198487 && (lumi >= 760 && lumi <= 852))
    return true;
  if (run == 198487 && (lumi >= 854 && lumi <= 994))
    return true;
  if (run == 198487 && (lumi >= 997 && lumi <= 1434))
    return true;
  if (run == 198487 && (lumi >= 1437 && lumi <= 1610))
    return true;
  if (run == 204554 && (lumi >= 1 && lumi <= 5))
    return true;
  if (run == 204554 && (lumi >= 7 && lumi <= 221))
    return true;
  if (run == 204554 && (lumi >= 224 && lumi <= 455))
    return true;
  if (run == 204554 && (lumi >= 458 && lumi <= 470))
    return true;
  if (run == 204554 && (lumi >= 472 && lumi <= 481))
    return true;
  if (run == 204554 && (lumi >= 483 && lumi <= 514))
    return true;
  if (run == 199336 && (lumi >= 1 && lumi <= 33))
    return true;
  if (run == 199336 && (lumi >= 36 && lumi <= 122))
    return true;
  if (run == 199336 && (lumi >= 125 && lumi <= 231))
    return true;
  if (run == 199336 && (lumi >= 234 && lumi <= 614))
    return true;
  if (run == 199336 && (lumi >= 617 && lumi <= 789))
    return true;
  if (run == 199336 && (lumi >= 791 && lumi <= 977))
    return true;
  if (run == 204552 && (lumi >= 1 && lumi <= 9))
    return true;
  if (run == 204553 && (lumi >= 1 && lumi <= 51))
    return true;
  if (run == 204553 && (lumi >= 53 && lumi <= 60))
    return true;
  if (run == 204553 && (lumi >= 63 && lumi <= 101))
    return true;
  if (run == 206940 && (lumi >= 1 && lumi <= 151))
    return true;
  if (run == 206940 && (lumi >= 153 && lumi <= 153))
    return true;
  if (run == 206940 && (lumi >= 155 && lumi <= 298))
    return true;
  if (run == 206940 && (lumi >= 301 && lumi <= 382))
    return true;
  if (run == 206940 && (lumi >= 384 && lumi <= 712))
    return true;
  if (run == 206940 && (lumi >= 715 && lumi <= 803))
    return true;
  if (run == 206940 && (lumi >= 805 && lumi <= 960))
    return true;
  if (run == 206940 && (lumi >= 963 && lumi <= 1027))
    return true;
  if (run == 193834 && (lumi >= 1 && lumi <= 35))
    return true;
  if (run == 193835 && (lumi >= 1 && lumi <= 20))
    return true;
  if (run == 193835 && (lumi >= 22 && lumi <= 26))
    return true;
  if (run == 193836 && (lumi >= 1 && lumi <= 2))
    return true;
  if (run == 202504 && (lumi >= 77 && lumi <= 96))
    return true;
  if (run == 202504 && (lumi >= 99 && lumi <= 133))
    return true;
  if (run == 202504 && (lumi >= 135 && lumi <= 182))
    return true;
  if (run == 202504 && (lumi >= 184 && lumi <= 211))
    return true;
  if (run == 202504 && (lumi >= 213 && lumi <= 241))
    return true;
  if (run == 202504 && (lumi >= 243 && lumi <= 392))
    return true;
  if (run == 202504 && (lumi >= 395 && lumi <= 527))
    return true;
  if (run == 202504 && (lumi >= 529 && lumi <= 617))
    return true;
  if (run == 202504 && (lumi >= 620 && lumi <= 715))
    return true;
  if (run == 202504 && (lumi >= 718 && lumi <= 763))
    return true;
  if (run == 202504 && (lumi >= 766 && lumi <= 1172))
    return true;
  if (run == 202504 && (lumi >= 1174 && lumi <= 1247))
    return true;
  if (run == 202504 && (lumi >= 1250 && lumi <= 1471))
    return true;
  if (run == 202504 && (lumi >= 1474 && lumi <= 1679))
    return true;
  if (run == 202504 && (lumi >= 1682 && lumi <= 1704))
    return true;
  if (run == 190646 && (lumi >= 1 && lumi <= 111))
    return true;
  if (run == 190645 && (lumi >= 10 && lumi <= 110))
    return true;
  if (run == 208427 && (lumi >= 49 && lumi <= 89))
    return true;
  if (run == 208427 && (lumi >= 92 && lumi <= 161))
    return true;
  if (run == 208427 && (lumi >= 164 && lumi <= 164))
    return true;
  if (run == 208427 && (lumi >= 166 && lumi <= 173))
    return true;
  if (run == 208427 && (lumi >= 175 && lumi <= 268))
    return true;
  if (run == 208427 && (lumi >= 271 && lumi <= 312))
    return true;
  if (run == 208427 && (lumi >= 315 && lumi <= 315))
    return true;
  if (run == 208427 && (lumi >= 317 && lumi <= 335))
    return true;
  if (run == 208427 && (lumi >= 337 && lumi <= 361))
    return true;
  if (run == 208427 && (lumi >= 364 && lumi <= 402))
    return true;
  if (run == 208427 && (lumi >= 404 && lumi <= 422))
    return true;
  if (run == 208427 && (lumi >= 425 && lumi <= 577))
    return true;
  if (run == 208427 && (lumi >= 580 && lumi <= 647))
    return true;
  if (run == 193999 && (lumi >= 1 && lumi <= 45))
    return true;
  if (run == 193998 && (lumi >= 66 && lumi <= 113))
    return true;
  if (run == 193998 && (lumi >= 115 && lumi <= 278))
    return true;
  if (run == 198969 && (lumi >= 58 && lumi <= 81))
    return true;
  if (run == 198969 && (lumi >= 84 && lumi <= 247))
    return true;
  if (run == 198969 && (lumi >= 249 && lumi <= 323))
    return true;
  if (run == 198969 && (lumi >= 325 && lumi <= 365))
    return true;
  if (run == 198969 && (lumi >= 367 && lumi <= 413))
    return true;
  if (run == 198969 && (lumi >= 416 && lumi <= 466))
    return true;
  if (run == 198969 && (lumi >= 468 && lumi <= 643))
    return true;
  if (run == 198969 && (lumi >= 646 && lumi <= 918))
    return true;
  if (run == 198969 && (lumi >= 920 && lumi <= 1011))
    return true;
  if (run == 198969 && (lumi >= 1013 && lumi <= 1175))
    return true;
  if (run == 198969 && (lumi >= 1178 && lumi <= 1236))
    return true;
  if (run == 198969 && (lumi >= 1239 && lumi <= 1253))
    return true;
  if (run == 208429 && (lumi >= 1 && lumi <= 56))
    return true;
  if (run == 208429 && (lumi >= 59 && lumi <= 139))
    return true;
  if (run == 208429 && (lumi >= 141 && lumi <= 159))
    return true;
  if (run == 208429 && (lumi >= 162 && lumi <= 237))
    return true;
  if (run == 208429 && (lumi >= 240 && lumi <= 440))
    return true;
  if (run == 208429 && (lumi >= 442 && lumi <= 452))
    return true;
  if (run == 208429 && (lumi >= 455 && lumi <= 589))
    return true;
  if (run == 208429 && (lumi >= 592 && lumi <= 712))
    return true;
  if (run == 208429 && (lumi >= 715 && lumi <= 922))
    return true;
  if (run == 206088 && (lumi >= 86 && lumi <= 159))
    return true;
  if (run == 206088 && (lumi >= 161 && lumi <= 178))
    return true;
  if (run == 206088 && (lumi >= 181 && lumi <= 199))
    return true;
  if (run == 206088 && (lumi >= 202 && lumi <= 286))
    return true;
  if (run == 195647 && (lumi >= 1 && lumi <= 41))
    return true;
  if (run == 201159 && (lumi >= 70 && lumi <= 211))
    return true;
  if (run == 195649 && (lumi >= 1 && lumi <= 69))
    return true;
  if (run == 195649 && (lumi >= 72 && lumi <= 151))
    return true;
  if (run == 195649 && (lumi >= 154 && lumi <= 181))
    return true;
  if (run == 195649 && (lumi >= 183 && lumi <= 247))
    return true;
  if (run == 200976 && (lumi >= 94 && lumi <= 164))
    return true;
  if (run == 202477 && (lumi >= 1 && lumi <= 129))
    return true;
  if (run == 202477 && (lumi >= 131 && lumi <= 150))
    return true;
  if (run == 202472 && (lumi >= 1 && lumi <= 96))
    return true;
  if (run == 202472 && (lumi >= 99 && lumi <= 112))
    return true;
  if (run == 204541 && (lumi >= 5 && lumi <= 39))
    return true;
  if (run == 204541 && (lumi >= 42 && lumi <= 42))
    return true;
  if (run == 204541 && (lumi >= 44 && lumi <= 139))
    return true;
  if (run == 204541 && (lumi >= 142 && lumi <= 149))
    return true;
  if (run == 204541 && (lumi >= 151 && lumi <= 204))
    return true;
  if (run == 202478 && (lumi >= 1 && lumi <= 177))
    return true;
  if (run == 202478 && (lumi >= 180 && lumi <= 183))
    return true;
  if (run == 202478 && (lumi >= 186 && lumi <= 219))
    return true;
  if (run == 202478 && (lumi >= 222 && lumi <= 360))
    return true;
  if (run == 202478 && (lumi >= 362 && lumi <= 506))
    return true;
  if (run == 202478 && (lumi >= 509 && lumi <= 531))
    return true;
  if (run == 202478 && (lumi >= 534 && lumi <= 718))
    return true;
  if (run == 202478 && (lumi >= 720 && lumi <= 927))
    return true;
  if (run == 202478 && (lumi >= 929 && lumi <= 973))
    return true;
  if (run == 202478 && (lumi >= 975 && lumi <= 1029))
    return true;
  if (run == 202478 && (lumi >= 1031 && lumi <= 1186))
    return true;
  if (run == 202478 && (lumi >= 1189 && lumi <= 1212))
    return true;
  if (run == 202478 && (lumi >= 1215 && lumi <= 1248))
    return true;
  if (run == 193193 && (lumi >= 1 && lumi <= 6))
    return true;
  if (run == 193193 && (lumi >= 8 && lumi <= 8))
    return true;
  if (run == 193193 && (lumi >= 11 && lumi <= 83))
    return true;
  if (run == 193193 && (lumi >= 86 && lumi <= 120))
    return true;
  if (run == 193193 && (lumi >= 122 && lumi <= 160))
    return true;
  if (run == 193193 && (lumi >= 162 && lumi <= 274))
    return true;
  if (run == 193193 && (lumi >= 276 && lumi <= 495))
    return true;
  if (run == 193193 && (lumi >= 497 && lumi <= 506))
    return true;
  if (run == 193192 && (lumi >= 58 && lumi <= 86))
    return true;
  if (run == 207397 && (lumi >= 32 && lumi <= 77))
    return true;
  if (run == 207397 && (lumi >= 80 && lumi <= 140))
    return true;
  if (run == 207397 && (lumi >= 143 && lumi <= 179))
    return true;
  if (run == 202093 && (lumi >= 1 && lumi <= 104))
    return true;
  if (run == 202093 && (lumi >= 107 && lumi <= 320))
    return true;
  if (run == 202093 && (lumi >= 322 && lumi <= 360))
    return true;
  if (run == 207477 && (lumi >= 76 && lumi <= 104))
    return true;
  if (run == 207477 && (lumi >= 107 && lumi <= 111))
    return true;
  if (run == 207477 && (lumi >= 114 && lumi <= 147))
    return true;
  if (run == 207477 && (lumi >= 150 && lumi <= 295))
    return true;
  if (run == 207477 && (lumi >= 298 && lumi <= 483))
    return true;
  if (run == 207477 && (lumi >= 486 && lumi <= 494))
    return true;
  if (run == 207477 && (lumi >= 497 && lumi <= 527))
    return true;
  if (run == 207477 && (lumi >= 530 && lumi <= 563))
    return true;
  if (run == 207477 && (lumi >= 565 && lumi <= 570))
    return true;
  if (run == 194424 && (lumi >= 63 && lumi <= 141))
    return true;
  if (run == 194424 && (lumi >= 144 && lumi <= 195))
    return true;
  if (run == 194424 && (lumi >= 198 && lumi <= 266))
    return true;
  if (run == 194424 && (lumi >= 268 && lumi <= 421))
    return true;
  if (run == 194424 && (lumi >= 424 && lumi <= 478))
    return true;
  if (run == 194424 && (lumi >= 481 && lumi <= 531))
    return true;
  if (run == 194424 && (lumi >= 534 && lumi <= 553))
    return true;
  if (run == 194424 && (lumi >= 556 && lumi <= 706))
    return true;
  if (run == 194424 && (lumi >= 708 && lumi <= 708))
    return true;
  if (run == 194428 && (lumi >= 1 && lumi <= 85))
    return true;
  if (run == 194428 && (lumi >= 87 && lumi <= 122))
    return true;
  if (run == 194428 && (lumi >= 125 && lumi <= 294))
    return true;
  if (run == 194428 && (lumi >= 296 && lumi <= 465))
    return true;
  if (run == 194429 && (lumi >= 1 && lumi <= 4))
    return true;
  if (run == 194429 && (lumi >= 7 && lumi <= 54))
    return true;
  if (run == 194429 && (lumi >= 57 && lumi <= 147))
    return true;
  if (run == 194429 && (lumi >= 150 && lumi <= 411))
    return true;
  if (run == 194429 && (lumi >= 413 && lumi <= 742))
    return true;
  if (run == 194429 && (lumi >= 745 && lumi <= 986))
    return true;
  if (run == 194429 && (lumi >= 988 && lumi <= 1019))
    return true;
  if (run == 199021 && (lumi >= 59 && lumi <= 88))
    return true;
  if (run == 199021 && (lumi >= 91 && lumi <= 128))
    return true;
  if (run == 199021 && (lumi >= 130 && lumi <= 133))
    return true;
  if (run == 199021 && (lumi >= 136 && lumi <= 309))
    return true;
  if (run == 199021 && (lumi >= 311 && lumi <= 333))
    return true;
  if (run == 199021 && (lumi >= 335 && lumi <= 410))
    return true;
  if (run == 199021 && (lumi >= 414 && lumi <= 469))
    return true;
  if (run == 199021 && (lumi >= 471 && lumi <= 533))
    return true;
  if (run == 199021 && (lumi >= 535 && lumi <= 563))
    return true;
  if (run == 199021 && (lumi >= 565 && lumi <= 1223))
    return true;
  if (run == 199021 && (lumi >= 1226 && lumi <= 1479))
    return true;
  if (run == 199021 && (lumi >= 1481 && lumi <= 1494))
    return true;
  if (run == 201278 && (lumi >= 62 && lumi <= 163))
    return true;
  if (run == 201278 && (lumi >= 166 && lumi <= 229))
    return true;
  if (run == 201278 && (lumi >= 232 && lumi <= 256))
    return true;
  if (run == 201278 && (lumi >= 259 && lumi <= 316))
    return true;
  if (run == 201278 && (lumi >= 318 && lumi <= 595))
    return true;
  if (run == 201278 && (lumi >= 598 && lumi <= 938))
    return true;
  if (run == 201278 && (lumi >= 942 && lumi <= 974))
    return true;
  if (run == 201278 && (lumi >= 976 && lumi <= 1160))
    return true;
  if (run == 201278 && (lumi >= 1163 && lumi <= 1304))
    return true;
  if (run == 201278 && (lumi >= 1306 && lumi <= 1793))
    return true;
  if (run == 201278 && (lumi >= 1796 && lumi <= 1802))
    return true;
  if (run == 201278 && (lumi >= 1805 && lumi <= 1906))
    return true;
  if (run == 201278 && (lumi >= 1909 && lumi <= 1929))
    return true;
  if (run == 201278 && (lumi >= 1932 && lumi <= 2174))
    return true;
  if (run == 205158 && (lumi >= 81 && lumi <= 289))
    return true;
  if (run == 205158 && (lumi >= 292 && lumi <= 313))
    return true;
  if (run == 205158 && (lumi >= 315 && lumi <= 473))
    return true;
  if (run == 205158 && (lumi >= 476 && lumi <= 591))
    return true;
  if (run == 205158 && (lumi >= 594 && lumi <= 595))
    return true;
  if (run == 205158 && (lumi >= 597 && lumi <= 612))
    return true;
  if (run == 205158 && (lumi >= 615 && lumi <= 663))
    return true;
  if (run == 205158 && (lumi >= 665 && lumi <= 667))
    return true;
  if (run == 205158 && (lumi >= 672 && lumi <= 685))
    return true;
  if (run == 205158 && (lumi >= 687 && lumi <= 733))
    return true;
  if (run == 191226 && (lumi >= 77 && lumi <= 78))
    return true;
  if (run == 191226 && (lumi >= 81 && lumi <= 831))
    return true;
  if (run == 191226 && (lumi >= 833 && lumi <= 1454))
    return true;
  if (run == 191226 && (lumi >= 1456 && lumi <= 1466))
    return true;
  if (run == 191226 && (lumi >= 1469 && lumi <= 1507))
    return true;
  if (run == 191226 && (lumi >= 1510 && lumi <= 1686))
    return true;
  if (run == 205774 && (lumi >= 1 && lumi <= 80))
    return true;
  if (run == 205777 && (lumi >= 1 && lumi <= 8))
    return true;
  if (run == 199812 && (lumi >= 70 && lumi <= 141))
    return true;
  if (run == 199812 && (lumi >= 144 && lumi <= 163))
    return true;
  if (run == 199812 && (lumi >= 182 && lumi <= 211))
    return true;
  if (run == 199812 && (lumi >= 214 && lumi <= 471))
    return true;
  if (run == 199812 && (lumi >= 474 && lumi <= 505))
    return true;
  if (run == 199812 && (lumi >= 508 && lumi <= 557))
    return true;
  if (run == 199812 && (lumi >= 560 && lumi <= 571))
    return true;
  if (run == 199812 && (lumi >= 574 && lumi <= 623))
    return true;
  if (run == 199812 && (lumi >= 626 && lumi <= 751))
    return true;
  if (run == 199812 && (lumi >= 754 && lumi <= 796))
    return true;
  if (run == 195868 && (lumi >= 1 && lumi <= 88))
    return true;
  if (run == 195868 && (lumi >= 90 && lumi <= 107))
    return true;
  if (run == 195868 && (lumi >= 110 && lumi <= 205))
    return true;
  if (run == 196363 && (lumi >= 1 && lumi <= 8))
    return true;
  if (run == 196363 && (lumi >= 11 && lumi <= 34))
    return true;
  if (run == 196362 && (lumi >= 1 && lumi <= 88))
    return true;
  if (run == 195399 && (lumi >= 1 && lumi <= 192))
    return true;
  if (run == 195399 && (lumi >= 194 && lumi <= 382))
    return true;
  if (run == 195398 && (lumi >= 3 && lumi <= 137))
    return true;
  if (run == 195398 && (lumi >= 139 && lumi <= 494))
    return true;
  if (run == 195398 && (lumi >= 497 && lumi <= 585))
    return true;
  if (run == 195398 && (lumi >= 587 && lumi <= 817))
    return true;
  if (run == 195398 && (lumi >= 820 && lumi <= 824))
    return true;
  if (run == 195398 && (lumi >= 827 && lumi <= 1225))
    return true;
  if (run == 195398 && (lumi >= 1228 && lumi <= 1307))
    return true;
  if (run == 195398 && (lumi >= 1309 && lumi <= 1712))
    return true;
  if (run == 195398 && (lumi >= 1721 && lumi <= 1736))
    return true;
  if (run == 195398 && (lumi >= 1741 && lumi <= 1752))
    return true;
  if (run == 195398 && (lumi >= 1767 && lumi <= 1795))
    return true;
  if (run == 201679 && (lumi >= 1 && lumi <= 110))
    return true;
  if (run == 201679 && (lumi >= 112 && lumi <= 241))
    return true;
  if (run == 201679 && (lumi >= 244 && lumi <= 298))
    return true;
  if (run == 201679 && (lumi >= 302 && lumi <= 321))
    return true;
  if (run == 201679 && (lumi >= 324 && lumi <= 461))
    return true;
  if (run == 201679 && (lumi >= 463 && lumi <= 483))
    return true;
  if (run == 194631 && (lumi >= 1 && lumi <= 42))
    return true;
  if (run == 194631 && (lumi >= 44 && lumi <= 100))
    return true;
  if (run == 194631 && (lumi >= 102 && lumi <= 169))
    return true;
  if (run == 194631 && (lumi >= 171 && lumi <= 222))
    return true;
  if (run == 194027 && (lumi >= 57 && lumi <= 113))
    return true;
  if (run == 195390 && (lumi >= 1 && lumi <= 1))
    return true;
  if (run == 195390 && (lumi >= 4 && lumi <= 27))
    return true;
  if (run == 195390 && (lumi >= 30 && lumi <= 145))
    return true;
  if (run == 195390 && (lumi >= 147 && lumi <= 183))
    return true;
  if (run == 195390 && (lumi >= 186 && lumi <= 187))
    return true;
  if (run == 195390 && (lumi >= 190 && lumi <= 208))
    return true;
  if (run == 195390 && (lumi >= 210 && lumi <= 213))
    return true;
  if (run == 195390 && (lumi >= 215 && lumi <= 400))
    return true;
  if (run == 195397 && (lumi >= 1 && lumi <= 10))
    return true;
  if (run == 195397 && (lumi >= 12 && lumi <= 89))
    return true;
  if (run == 195397 && (lumi >= 92 && lumi <= 120))
    return true;
  if (run == 195397 && (lumi >= 123 && lumi <= 141))
    return true;
  if (run == 195397 && (lumi >= 143 && lumi <= 251))
    return true;
  if (run == 195397 && (lumi >= 253 && lumi <= 253))
    return true;
  if (run == 195397 && (lumi >= 256 && lumi <= 475))
    return true;
  if (run == 195397 && (lumi >= 478 && lumi <= 525))
    return true;
  if (run == 195397 && (lumi >= 527 && lumi <= 608))
    return true;
  if (run == 195397 && (lumi >= 611 && lumi <= 776))
    return true;
  if (run == 195397 && (lumi >= 779 && lumi <= 970))
    return true;
  if (run == 195397 && (lumi >= 972 && lumi <= 1121))
    return true;
  if (run == 195397 && (lumi >= 1123 && lumi <= 1181))
    return true;
  if (run == 195397 && (lumi >= 1184 && lumi <= 1198))
    return true;
  if (run == 195397 && (lumi >= 1200 && lumi <= 1209))
    return true;
  if (run == 195396 && (lumi >= 49 && lumi <= 55))
    return true;
  if (run == 195396 && (lumi >= 58 && lumi <= 63))
    return true;
  if (run == 195396 && (lumi >= 66 && lumi <= 131))
    return true;
  if (run == 194735 && (lumi >= 44 && lumi <= 71))
    return true;
  if (run == 194735 && (lumi >= 74 && lumi <= 101))
    return true;
  if (run == 194735 && (lumi >= 104 && lumi <= 130))
    return true;
  if (run == 195551 && (lumi >= 91 && lumi <= 106))
    return true;
  if (run == 195552 && (lumi >= 1 && lumi <= 21))
    return true;
  if (run == 195552 && (lumi >= 23 && lumi <= 27))
    return true;
  if (run == 195552 && (lumi >= 30 && lumi <= 147))
    return true;
  if (run == 195552 && (lumi >= 149 && lumi <= 155))
    return true;
  if (run == 195552 && (lumi >= 158 && lumi <= 182))
    return true;
  if (run == 195552 && (lumi >= 185 && lumi <= 287))
    return true;
  if (run == 195552 && (lumi >= 290 && lumi <= 349))
    return true;
  if (run == 195552 && (lumi >= 352 && lumi <= 469))
    return true;
  if (run == 195552 && (lumi >= 472 && lumi <= 815))
    return true;
  if (run == 195552 && (lumi >= 818 && lumi <= 823))
    return true;
  if (run == 195552 && (lumi >= 825 && lumi <= 883))
    return true;
  if (run == 195552 && (lumi >= 885 && lumi <= 1152))
    return true;
  if (run == 195552 && (lumi >= 1154 && lumi <= 1300))
    return true;
  if (run == 195552 && (lumi >= 1303 && lumi <= 1789))
    return true;
  if (run == 203002 && (lumi >= 77 && lumi <= 128))
    return true;
  if (run == 203002 && (lumi >= 130 && lumi <= 141))
    return true;
  if (run == 203002 && (lumi >= 144 && lumi <= 207))
    return true;
  if (run == 203002 && (lumi >= 209 && lumi <= 267))
    return true;
  if (run == 203002 && (lumi >= 270 && lumi <= 360))
    return true;
  if (run == 203002 && (lumi >= 362 && lumi <= 501))
    return true;
  if (run == 203002 && (lumi >= 504 && lumi <= 641))
    return true;
  if (run == 203002 && (lumi >= 643 && lumi <= 669))
    return true;
  if (run == 203002 && (lumi >= 671 && lumi <= 671))
    return true;
  if (run == 203002 && (lumi >= 674 && lumi <= 717))
    return true;
  if (run == 203002 && (lumi >= 720 && lumi <= 1034))
    return true;
  if (run == 203002 && (lumi >= 1037 && lumi <= 1070))
    return true;
  if (run == 203002 && (lumi >= 1073 && lumi <= 1370))
    return true;
  if (run == 203002 && (lumi >= 1372 && lumi <= 1392))
    return true;
  if (run == 203002 && (lumi >= 1395 && lumi <= 1410))
    return true;
  if (run == 203002 && (lumi >= 1413 && lumi <= 1596))
    return true;
  if (run == 201613 && (lumi >= 1 && lumi <= 42))
    return true;
  if (run == 201613 && (lumi >= 44 && lumi <= 49))
    return true;
  if (run == 201613 && (lumi >= 53 && lumi <= 210))
    return true;
  if (run == 201613 && (lumi >= 213 && lumi <= 215))
    return true;
  if (run == 201613 && (lumi >= 218 && lumi <= 225))
    return true;
  if (run == 201613 && (lumi >= 228 && lumi <= 646))
    return true;
  if (run == 201611 && (lumi >= 87 && lumi <= 145))
    return true;
  if (run == 201611 && (lumi >= 149 && lumi <= 182))
    return true;
  if (run == 201611 && (lumi >= 184 && lumi <= 186))
    return true;
  if (run == 195919 && (lumi >= 1 && lumi <= 15))
    return true;
  if (run == 195918 && (lumi >= 1 && lumi <= 44))
    return true;
  if (run == 195918 && (lumi >= 46 && lumi <= 46))
    return true;
  if (run == 195918 && (lumi >= 49 && lumi <= 64))
    return true;
  if (run == 195915 && (lumi >= 1 && lumi <= 109))
    return true;
  if (run == 195915 && (lumi >= 111 && lumi <= 275))
    return true;
  if (run == 195915 && (lumi >= 278 && lumi <= 390))
    return true;
  if (run == 195915 && (lumi >= 393 && lumi <= 417))
    return true;
  if (run == 195915 && (lumi >= 419 && lumi <= 429))
    return true;
  if (run == 195915 && (lumi >= 432 && lumi <= 505))
    return true;
  if (run == 195915 && (lumi >= 507 && lumi <= 747))
    return true;
  if (run == 195915 && (lumi >= 749 && lumi <= 785))
    return true;
  if (run == 195915 && (lumi >= 787 && lumi <= 828))
    return true;
  if (run == 195915 && (lumi >= 830 && lumi <= 850))
    return true;
  if (run == 195917 && (lumi >= 1 && lumi <= 4))
    return true;
  if (run == 195916 && (lumi >= 1 && lumi <= 16))
    return true;
  if (run == 195916 && (lumi >= 19 && lumi <= 68))
    return true;
  if (run == 195916 && (lumi >= 71 && lumi <= 212))
    return true;
  if (run == 200599 && (lumi >= 75 && lumi <= 129))
    return true;
  if (run == 200599 && (lumi >= 132 && lumi <= 137))
    return true;
  if (run == 206572 && (lumi >= 37 && lumi <= 47))
    return true;
  if (run == 206573 && (lumi >= 2 && lumi <= 14))
    return true;
  if (run == 206478 && (lumi >= 1 && lumi <= 27))
    return true;
  if (run == 206478 && (lumi >= 29 && lumi <= 136))
    return true;
  if (run == 206478 && (lumi >= 139 && lumi <= 144))
    return true;
  if (run == 206574 && (lumi >= 1 && lumi <= 87))
    return true;
  if (run == 206575 && (lumi >= 1 && lumi <= 7))
    return true;
  if (run == 206575 && (lumi >= 10 && lumi <= 10))
    return true;
  if (run == 206575 && (lumi >= 12 && lumi <= 69))
    return true;
  if (run == 206477 && (lumi >= 1 && lumi <= 14))
    return true;
  if (run == 206477 && (lumi >= 16 && lumi <= 31))
    return true;
  if (run == 206477 && (lumi >= 33 && lumi <= 41))
    return true;
  if (run == 206477 && (lumi >= 44 && lumi <= 51))
    return true;
  if (run == 206477 && (lumi >= 53 && lumi <= 70))
    return true;
  if (run == 206477 && (lumi >= 73 && lumi <= 75))
    return true;
  if (run == 206477 && (lumi >= 77 && lumi <= 89))
    return true;
  if (run == 206477 && (lumi >= 91 && lumi <= 94))
    return true;
  if (run == 206477 && (lumi >= 97 && lumi <= 115))
    return true;
  if (run == 206477 && (lumi >= 118 && lumi <= 184))
    return true;
  if (run == 206476 && (lumi >= 73 && lumi <= 129))
    return true;
  if (run == 206476 && (lumi >= 133 && lumi <= 137))
    return true;
  if (run == 206476 && (lumi >= 140 && lumi <= 141))
    return true;
  if (run == 206476 && (lumi >= 143 && lumi <= 219))
    return true;
  if (run == 207372 && (lumi >= 1 && lumi <= 27))
    return true;
  if (run == 207372 && (lumi >= 30 && lumi <= 113))
    return true;
  if (run == 207372 && (lumi >= 116 && lumi <= 154))
    return true;
  if (run == 207372 && (lumi >= 156 && lumi <= 174))
    return true;
  if (run == 207372 && (lumi >= 176 && lumi <= 478))
    return true;
  if (run == 207372 && (lumi >= 480 && lumi <= 496))
    return true;
  if (run == 202016 && (lumi >= 1 && lumi <= 48))
    return true;
  if (run == 202016 && (lumi >= 51 && lumi <= 134))
    return true;
  if (run == 202016 && (lumi >= 137 && lumi <= 177))
    return true;
  if (run == 202016 && (lumi >= 179 && lumi <= 743))
    return true;
  if (run == 202016 && (lumi >= 745 && lumi <= 831))
    return true;
  if (run == 202016 && (lumi >= 834 && lumi <= 890))
    return true;
  if (run == 202016 && (lumi >= 893 && lumi <= 896))
    return true;
  if (run == 202016 && (lumi >= 898 && lumi <= 932))
    return true;
  if (run == 202016 && (lumi >= 934 && lumi <= 1010))
    return true;
  if (run == 207273 && (lumi >= 3 && lumi <= 877))
    return true;
  if (run == 190782 && (lumi >= 55 && lumi <= 181))
    return true;
  if (run == 190782 && (lumi >= 184 && lumi <= 233))
    return true;
  if (run == 190782 && (lumi >= 236 && lumi <= 399))
    return true;
  if (run == 190782 && (lumi >= 401 && lumi <= 409))
    return true;
  if (run == 206901 && (lumi >= 1 && lumi <= 98))
    return true;
  if (run == 204511 && (lumi >= 1 && lumi <= 56))
    return true;
  if (run == 206906 && (lumi >= 1 && lumi <= 31))
    return true;
  if (run == 206906 && (lumi >= 38 && lumi <= 94))
    return true;
  if (run == 206906 && (lumi >= 96 && lumi <= 136))
    return true;
  if (run == 206906 && (lumi >= 138 && lumi <= 139))
    return true;
  if (run == 206906 && (lumi >= 142 && lumi <= 149))
    return true;
  if (run == 206906 && (lumi >= 151 && lumi <= 175))
    return true;
  if (run == 206906 && (lumi >= 177 && lumi <= 206))
    return true;
  if (run == 207371 && (lumi >= 72 && lumi <= 117))
    return true;
  if (run == 207371 && (lumi >= 120 && lumi <= 124))
    return true;
  if (run == 207279 && (lumi >= 68 && lumi <= 138))
    return true;
  if (run == 207279 && (lumi >= 141 && lumi <= 149))
    return true;
  if (run == 207279 && (lumi >= 151 && lumi <= 237))
    return true;
  if (run == 207279 && (lumi >= 240 && lumi <= 266))
    return true;
  if (run == 207279 && (lumi >= 269 && lumi <= 307))
    return true;
  if (run == 207279 && (lumi >= 309 && lumi <= 416))
    return true;
  if (run == 207279 && (lumi >= 498 && lumi <= 551))
    return true;
  if (run == 207279 && (lumi >= 554 && lumi <= 640))
    return true;
  if (run == 207279 && (lumi >= 643 && lumi <= 961))
    return true;
  if (run == 207279 && (lumi >= 963 && lumi <= 1095))
    return true;
  if (run == 207279 && (lumi >= 1098 && lumi <= 1160))
    return true;
  return false;
}

} // namespace LumiFilter
