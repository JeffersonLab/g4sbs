/*
 *  cteq-database.c
 *  
 *
 *  Created by Zoltan Nagy on 3/22/08.
 *  Copyright 2008 Zoltan Nagy. All rights reserved.
 *
 */

/* Standard C includes */
#include <string.h>


/* cteq-pdf includes */
#include "cteqpdf.h"

#define CTEQ6STD_TBL_PATH    CTEQ_TBL_PATH"/cteq6std/"
#define CTEQ6_TBL_PATH       CTEQ_TBL_PATH"/cteq6/"
#define CTEQ65S_TBL_PATH     CTEQ_TBL_PATH"/ctq65s/"
#define CTEQ65C_TBL_PATH     CTEQ_TBL_PATH"/ctq65c/"
#define CTEQ6M_TBL_PATH      CTEQ_TBL_PATH"/cteq6m/"
#define CTEQ61_TBL_PATH      CTEQ_TBL_PATH"/cteq61/"
#define CTEQ65_PDS_TBL_PATH  CTEQ_TBL_PATH"/ctq65-pds/"
#define CTEQ66A_TBL_PATH     CTEQ_TBL_PATH"/ctq66a/"
#define CTEQ66C_TBL_PATH     CTEQ_TBL_PATH"/ctq66c/"
#define CTEQ66M_TBL_PATH     CTEQ_TBL_PATH"/ctq66m/"



static const cteq_pdfset_t __cteq_pdfset_database[] = { 
  {1,   "CTEQ6M",     "description", CTEQ6STD_TBL_PATH, "cteq6m.tbl",  1},
  {2,   "CTEQ6D",     "description", CTEQ6STD_TBL_PATH, "cteq6d.tbl",  1},
  {3,   "CTEQ6L",     "description", CTEQ6STD_TBL_PATH, "cteq6l.tbl",  1},
  {4,   "CTEQ6L1",    "description", CTEQ6STD_TBL_PATH, "cteq6l1.tbl", 1},

  {11,  "CTEQ6A",     "description", CTEQ6_TBL_PATH, "cteq6sa.pds",  2},
  {12,  "CTEQ6B",     "description", CTEQ6_TBL_PATH, "cteq6sb.pds",  2},
  {13,  "CTEQ6C",     "description", CTEQ6_TBL_PATH, "cteq6sc.pds",  2},
  {14,  "CTEQ6B+",    "description", CTEQ6_TBL_PATH, "cteq6sb+.pds", 2},
  {15,  "CTEQ6B-",    "description", CTEQ6_TBL_PATH, "cteq6sb-.pds", 2},
  {21,  "CTEQ6HQ",    "description", CTEQ6_TBL_PATH, "cteq6hq.pds",  2},
 
  {30,  "CTEQ6.5S0",  "description", CTEQ65S_TBL_PATH, "ctq65.s+0.pds", 2},
  {31,  "CTEQ6.5S1",  "description", CTEQ65S_TBL_PATH, "ctq65.s+1.pds", 2},
  {32,  "CTEQ6.5S2",  "description", CTEQ65S_TBL_PATH, "ctq65.s+2.pds", 2},
  {33,  "CTEQ6.5S3",  "description", CTEQ65S_TBL_PATH, "ctq65.s+3.pds", 2},
  {34,  "CTEQ6.5S4",  "description", CTEQ65S_TBL_PATH, "ctq65.s+4.pds", 2},
  {35,  "CTEQ6.5S-0", "description", CTEQ65S_TBL_PATH, "ctq65.s-0.pds", 2},
  {36,  "CTEQ6.5S-1", "description", CTEQ65S_TBL_PATH, "ctq65.s-1.pds", 2},
  {37,  "CTEQ6.5S-2", "description", CTEQ65S_TBL_PATH, "ctq65.s-2.pds", 2},

  {40,  "CTEQ6.5C0",  "description", CTEQ65C_TBL_PATH, "ctq65.c0.pds", 2},
  {41,  "CTEQ6.5C1",  "description", CTEQ65C_TBL_PATH, "ctq65.c1.pds", 2},
  {42,  "CTEQ6.5C2",  "description", CTEQ65C_TBL_PATH, "ctq65.c2.pds", 2},
  {43,  "CTEQ6.5C3",  "description", CTEQ65C_TBL_PATH, "ctq65.c3.pds", 2},
  {44,  "CTEQ6.5C4",  "description", CTEQ65C_TBL_PATH, "ctq65.c4.pds", 2},
  {45,  "CTEQ6.5C5",  "description", CTEQ65C_TBL_PATH, "ctq65.c5.pds", 2},
  {46,  "CTEQ6.5C6",  "description", CTEQ65C_TBL_PATH, "ctq65.c6.pds", 2},

  {100, "CTEQ6M.00",  "description", CTEQ6M_TBL_PATH, "cteq6m100.tbl", 1},
  {101, "CTEQ6M.01",  "description", CTEQ6M_TBL_PATH, "cteq6m101.tbl", 1},
  {102, "CTEQ6M.02",  "description", CTEQ6M_TBL_PATH, "cteq6m102.tbl", 1},
  {103, "CTEQ6M.03",  "description", CTEQ6M_TBL_PATH, "cteq6m103.tbl", 1},
  {104, "CTEQ6M.04",  "description", CTEQ6M_TBL_PATH, "cteq6m104.tbl", 1},
  {105, "CTEQ6M.05",  "description", CTEQ6M_TBL_PATH, "cteq6m105.tbl", 1},
  {106, "CTEQ6M.06",  "description", CTEQ6M_TBL_PATH, "cteq6m106.tbl", 1},
  {107, "CTEQ6M.07",  "description", CTEQ6M_TBL_PATH, "cteq6m107.tbl", 1},
  {108, "CTEQ6M.08",  "description", CTEQ6M_TBL_PATH, "cteq6m108.tbl", 1},
  {109, "CTEQ6M.09",  "description", CTEQ6M_TBL_PATH, "cteq6m109.tbl", 1},
  {110, "CTEQ6M.10",  "description", CTEQ6M_TBL_PATH, "cteq6m110.tbl", 1},
  {111, "CTEQ6M.11",  "description", CTEQ6M_TBL_PATH, "cteq6m111.tbl", 1},
  {112, "CTEQ6M.12",  "description", CTEQ6M_TBL_PATH, "cteq6m112.tbl", 1},
  {113, "CTEQ6M.13",  "description", CTEQ6M_TBL_PATH, "cteq6m113.tbl", 1},
  {114, "CTEQ6M.14",  "description", CTEQ6M_TBL_PATH, "cteq6m114.tbl", 1},
  {115, "CTEQ6M.15",  "description", CTEQ6M_TBL_PATH, "cteq6m115.tbl", 1},
  {116, "CTEQ6M.16",  "description", CTEQ6M_TBL_PATH, "cteq6m116.tbl", 1},
  {117, "CTEQ6M.17",  "description", CTEQ6M_TBL_PATH, "cteq6m117.tbl", 1},
  {118, "CTEQ6M.18",  "description", CTEQ6M_TBL_PATH, "cteq6m118.tbl", 1},
  {119, "CTEQ6M.19",  "description", CTEQ6M_TBL_PATH, "cteq6m119.tbl", 1},
  {120, "CTEQ6M.20",  "description", CTEQ6M_TBL_PATH, "cteq6m120.tbl", 1},  
  {121, "CTEQ6M.21",  "description", CTEQ6M_TBL_PATH, "cteq6m121.tbl", 1},
  {122, "CTEQ6M.22",  "description", CTEQ6M_TBL_PATH, "cteq6m122.tbl", 1},
  {123, "CTEQ6M.23",  "description", CTEQ6M_TBL_PATH, "cteq6m123.tbl", 1},
  {124, "CTEQ6M.24",  "description", CTEQ6M_TBL_PATH, "cteq6m124.tbl", 1},
  {125, "CTEQ6M.25",  "description", CTEQ6M_TBL_PATH, "cteq6m125.tbl", 1},
  {126, "CTEQ6M.26",  "description", CTEQ6M_TBL_PATH, "cteq6m126.tbl", 1},
  {127, "CTEQ6M.27",  "description", CTEQ6M_TBL_PATH, "cteq6m127.tbl", 1},
  {128, "CTEQ6M.28",  "description", CTEQ6M_TBL_PATH, "cteq6m128.tbl", 1},
  {129, "CTEQ6M.29",  "description", CTEQ6M_TBL_PATH, "cteq6m129.tbl", 1},
  {130, "CTEQ6M.30",  "description", CTEQ6M_TBL_PATH, "cteq6m130.tbl", 1},
  {131, "CTEQ6M.31",  "description", CTEQ6M_TBL_PATH, "cteq6m131.tbl", 1},
  {132, "CTEQ6M.32",  "description", CTEQ6M_TBL_PATH, "cteq6m132.tbl", 1},
  {133, "CTEQ6M.33",  "description", CTEQ6M_TBL_PATH, "cteq6m133.tbl", 1},
  {134, "CTEQ6M.34",  "description", CTEQ6M_TBL_PATH, "cteq6m134.tbl", 1},
  {135, "CTEQ6M.35",  "description", CTEQ6M_TBL_PATH, "cteq6m135.tbl", 1},
  {136, "CTEQ6M.36",  "description", CTEQ6M_TBL_PATH, "cteq6m136.tbl", 1},
  {137, "CTEQ6M.37",  "description", CTEQ6M_TBL_PATH, "cteq6m137.tbl", 1},
  {138, "CTEQ6M.38",  "description", CTEQ6M_TBL_PATH, "cteq6m138.tbl", 1},
  {139, "CTEQ6M.39",  "description", CTEQ6M_TBL_PATH, "cteq6m139.tbl", 1},
  {140, "CTEQ6M.40",  "description", CTEQ6M_TBL_PATH, "cteq6m140.tbl", 1},
  
  {200, "CTEQ61.00",  "description", CTEQ61_TBL_PATH, "ctq61.00.tbl", 1},
  {201, "CTEQ61.01",  "description", CTEQ61_TBL_PATH, "ctq61.01.tbl", 1},
  {202, "CTEQ61.02",  "description", CTEQ61_TBL_PATH, "ctq61.02.tbl", 1},
  {203, "CTEQ61.03",  "description", CTEQ61_TBL_PATH, "ctq61.03.tbl", 1},
  {204, "CTEQ61.04",  "description", CTEQ61_TBL_PATH, "ctq61.04.tbl", 1},
  {205, "CTEQ61.05",  "description", CTEQ61_TBL_PATH, "ctq61.05.tbl", 1},
  {206, "CTEQ61.06",  "description", CTEQ61_TBL_PATH, "ctq61.06.tbl", 1},
  {207, "CTEQ61.07",  "description", CTEQ61_TBL_PATH, "ctq61.07.tbl", 1},
  {208, "CTEQ61.08",  "description", CTEQ61_TBL_PATH, "ctq61.08.tbl", 1},
  {209, "CTEQ61.09",  "description", CTEQ61_TBL_PATH, "ctq61.09.tbl", 1},
  {210, "CTEQ61.10",  "description", CTEQ61_TBL_PATH, "ctq61.10.tbl", 1},
  {211, "CTEQ61.11",  "description", CTEQ61_TBL_PATH, "ctq61.11.tbl", 1},
  {212, "CTEQ61.12",  "description", CTEQ61_TBL_PATH, "ctq61.12.tbl", 1},
  {213, "CTEQ61.13",  "description", CTEQ61_TBL_PATH, "ctq61.13.tbl", 1},
  {214, "CTEQ61.14",  "description", CTEQ61_TBL_PATH, "ctq61.14.tbl", 1},
  {215, "CTEQ61.15",  "description", CTEQ61_TBL_PATH, "ctq61.15.tbl", 1},
  {216, "CTEQ61.16",  "description", CTEQ61_TBL_PATH, "ctq61.16.tbl", 1},
  {217, "CTEQ61.17",  "description", CTEQ61_TBL_PATH, "ctq61.17.tbl", 1},
  {218, "CTEQ61.18",  "description", CTEQ61_TBL_PATH, "ctq61.18.tbl", 1},
  {219, "CTEQ61.19",  "description", CTEQ61_TBL_PATH, "ctq61.19.tbl", 1},
  {220, "CTEQ61.20",  "description", CTEQ61_TBL_PATH, "ctq61.20.tbl", 1},  
  {221, "CTEQ61.21",  "description", CTEQ61_TBL_PATH, "ctq61.21.tbl", 1},
  {222, "CTEQ61.22",  "description", CTEQ61_TBL_PATH, "ctq61.22.tbl", 1},
  {223, "CTEQ61.23",  "description", CTEQ61_TBL_PATH, "ctq61.23.tbl", 1},
  {224, "CTEQ61.24",  "description", CTEQ61_TBL_PATH, "ctq61.24.tbl", 1},
  {225, "CTEQ61.25",  "description", CTEQ61_TBL_PATH, "ctq61.25.tbl", 1},
  {226, "CTEQ61.26",  "description", CTEQ61_TBL_PATH, "ctq61.26.tbl", 1},
  {227, "CTEQ61.27",  "description", CTEQ61_TBL_PATH, "ctq61.27.tbl", 1},
  {228, "CTEQ61.28",  "description", CTEQ61_TBL_PATH, "ctq61.28.tbl", 1},
  {229, "CTEQ61.29",  "description", CTEQ61_TBL_PATH, "ctq61.29.tbl", 1},
  {230, "CTEQ61.30",  "description", CTEQ61_TBL_PATH, "ctq61.30.tbl", 1},
  {231, "CTEQ61.31",  "description", CTEQ61_TBL_PATH, "ctq61.31.tbl", 1},
  {232, "CTEQ61.32",  "description", CTEQ61_TBL_PATH, "ctq61.32.tbl", 1},
  {233, "CTEQ61.33",  "description", CTEQ61_TBL_PATH, "ctq61.33.tbl", 1},
  {234, "CTEQ61.34",  "description", CTEQ61_TBL_PATH, "ctq61.34.tbl", 1},
  {235, "CTEQ61.35",  "description", CTEQ61_TBL_PATH, "ctq61.35.tbl", 1},
  {236, "CTEQ61.36",  "description", CTEQ61_TBL_PATH, "ctq61.36.tbl", 1},
  {237, "CTEQ61.37",  "description", CTEQ61_TBL_PATH, "ctq61.37.tbl", 1},
  {238, "CTEQ61.38",  "description", CTEQ61_TBL_PATH, "ctq61.38.tbl", 1},
  {239, "CTEQ61.39",  "description", CTEQ61_TBL_PATH, "ctq61.39.tbl", 1},
  {240, "CTEQ61.40",  "description", CTEQ61_TBL_PATH, "ctq61.40.tbl", 1},
  
  {300, "CTEQ65.00",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.00.pds", 2},
  {301, "CTEQ65.01",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.01.pds", 2},
  {302, "CTEQ65.02",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.02.pds", 2},
  {303, "CTEQ65.03",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.03.pds", 2},
  {304, "CTEQ65.04",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.04.pds", 2},
  {305, "CTEQ65.05",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.05.pds", 2},
  {306, "CTEQ65.06",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.06.pds", 2},
  {307, "CTEQ65.07",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.07.pds", 2},
  {308, "CTEQ65.08",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.08.pds", 2},
  {309, "CTEQ65.09",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.09.pds", 2},
  {310, "CTEQ65.10",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.10.pds", 2},
  {311, "CTEQ65.11",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.11.pds", 2},
  {312, "CTEQ65.12",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.12.pds", 2},
  {313, "CTEQ65.13",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.13.pds", 2},
  {314, "CTEQ65.14",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.14.pds", 2},
  {315, "CTEQ65.15",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.15.pds", 2},
  {316, "CTEQ65.16",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.16.pds", 2},
  {317, "CTEQ65.17",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.17.pds", 2},
  {318, "CTEQ65.18",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.18.pds", 2},
  {319, "CTEQ65.19",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.19.pds", 2},
  {320, "CTEQ65.20",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.20.pds", 2},  
  {321, "CTEQ65.21",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.21.pds", 2},
  {322, "CTEQ65.22",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.22.pds", 2},
  {323, "CTEQ65.23",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.23.pds", 2},
  {324, "CTEQ65.24",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.24.pds", 2},
  {325, "CTEQ65.25",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.25.pds", 2},
  {326, "CTEQ65.26",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.26.pds", 2},
  {327, "CTEQ65.27",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.27.pds", 2},
  {328, "CTEQ65.28",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.28.pds", 2},
  {329, "CTEQ65.29",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.29.pds", 2},
  {330, "CTEQ65.30",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.30.pds", 2},
  {331, "CTEQ65.31",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.31.pds", 2},
  {332, "CTEQ65.32",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.32.pds", 2},
  {333, "CTEQ65.33",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.33.pds", 2},
  {334, "CTEQ65.34",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.34.pds", 2},
  {335, "CTEQ65.35",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.35.pds", 2},
  {336, "CTEQ65.36",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.36.pds", 2},
  {337, "CTEQ65.37",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.37.pds", 2},
  {338, "CTEQ65.38",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.38.pds", 2},
  {339, "CTEQ65.39",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.39.pds", 2},
  {340, "CTEQ65.40",  "description", CTEQ65_PDS_TBL_PATH, "ctq65.40.pds", 2},
  
  {400, "CTEQ66.00",  "description", CTEQ66M_TBL_PATH, "ctq66.00.pds", 2},
  {401, "CTEQ66.01",  "description", CTEQ66M_TBL_PATH, "ctq66.01.pds", 2},
  {402, "CTEQ66.02",  "description", CTEQ66M_TBL_PATH, "ctq66.02.pds", 2},
  {403, "CTEQ66.03",  "description", CTEQ66M_TBL_PATH, "ctq66.03.pds", 2},
  {404, "CTEQ66.04",  "description", CTEQ66M_TBL_PATH, "ctq66.04.pds", 2},
  {405, "CTEQ66.05",  "description", CTEQ66M_TBL_PATH, "ctq66.05.pds", 2},
  {406, "CTEQ66.06",  "description", CTEQ66M_TBL_PATH, "ctq66.06.pds", 2},
  {407, "CTEQ66.07",  "description", CTEQ66M_TBL_PATH, "ctq66.07.pds", 2},
  {408, "CTEQ66.08",  "description", CTEQ66M_TBL_PATH, "ctq66.08.pds", 2},
  {409, "CTEQ66.09",  "description", CTEQ66M_TBL_PATH, "ctq66.09.pds", 2},
  {410, "CTEQ66.10",  "description", CTEQ66M_TBL_PATH, "ctq66.10.pds", 2},
  {411, "CTEQ66.11",  "description", CTEQ66M_TBL_PATH, "ctq66.11.pds", 2},
  {412, "CTEQ66.12",  "description", CTEQ66M_TBL_PATH, "ctq66.12.pds", 2},
  {413, "CTEQ66.13",  "description", CTEQ66M_TBL_PATH, "ctq66.13.pds", 2},
  {414, "CTEQ66.14",  "description", CTEQ66M_TBL_PATH, "ctq66.14.pds", 2},
  {415, "CTEQ66.15",  "description", CTEQ66M_TBL_PATH, "ctq66.15.pds", 2},
  {416, "CTEQ66.16",  "description", CTEQ66M_TBL_PATH, "ctq66.16.pds", 2},
  {417, "CTEQ66.17",  "description", CTEQ66M_TBL_PATH, "ctq66.17.pds", 2},
  {418, "CTEQ66.18",  "description", CTEQ66M_TBL_PATH, "ctq66.18.pds", 2},
  {419, "CTEQ66.19",  "description", CTEQ66M_TBL_PATH, "ctq66.19.pds", 2},
  {420, "CTEQ66.20",  "description", CTEQ66M_TBL_PATH, "ctq66.20.pds", 2},  
  {421, "CTEQ66.21",  "description", CTEQ66M_TBL_PATH, "ctq66.21.pds", 2},
  {422, "CTEQ66.22",  "description", CTEQ66M_TBL_PATH, "ctq66.22.pds", 2},
  {423, "CTEQ66.23",  "description", CTEQ66M_TBL_PATH, "ctq66.23.pds", 2},
  {424, "CTEQ66.24",  "description", CTEQ66M_TBL_PATH, "ctq66.24.pds", 2},
  {425, "CTEQ66.25",  "description", CTEQ66M_TBL_PATH, "ctq66.25.pds", 2},
  {426, "CTEQ66.26",  "description", CTEQ66M_TBL_PATH, "ctq66.26.pds", 2},
  {427, "CTEQ66.27",  "description", CTEQ66M_TBL_PATH, "ctq66.27.pds", 2},
  {428, "CTEQ66.28",  "description", CTEQ66M_TBL_PATH, "ctq66.28.pds", 2},
  {429, "CTEQ66.29",  "description", CTEQ66M_TBL_PATH, "ctq66.29.pds", 2},
  {430, "CTEQ66.30",  "description", CTEQ66M_TBL_PATH, "ctq66.30.pds", 2},
  {431, "CTEQ66.31",  "description", CTEQ66M_TBL_PATH, "ctq66.31.pds", 2},
  {432, "CTEQ66.32",  "description", CTEQ66M_TBL_PATH, "ctq66.32.pds", 2},
  {433, "CTEQ66.33",  "description", CTEQ66M_TBL_PATH, "ctq66.33.pds", 2},
  {434, "CTEQ66.34",  "description", CTEQ66M_TBL_PATH, "ctq66.34.pds", 2},
  {435, "CTEQ66.35",  "description", CTEQ66M_TBL_PATH, "ctq66.35.pds", 2},
  {436, "CTEQ66.36",  "description", CTEQ66M_TBL_PATH, "ctq66.36.pds", 2},
  {437, "CTEQ66.37",  "description", CTEQ66M_TBL_PATH, "ctq66.37.pds", 2},
  {438, "CTEQ66.38",  "description", CTEQ66M_TBL_PATH, "ctq66.38.pds", 2},
  {439, "CTEQ66.39",  "description", CTEQ66M_TBL_PATH, "ctq66.39.pds", 2},
  {440, "CTEQ66.40",  "description", CTEQ66M_TBL_PATH, "ctq66.40.pds", 2},
  {441, "CTEQ66.41",  "description", CTEQ66M_TBL_PATH, "ctq66.41.pds", 2},
  {442, "CTEQ66.42",  "description", CTEQ66M_TBL_PATH, "ctq66.42.pds", 2},
  {443, "CTEQ66.43",  "description", CTEQ66M_TBL_PATH, "ctq66.43.pds", 2},
  {444, "CTEQ66.44",  "description", CTEQ66M_TBL_PATH, "ctq66.44.pds", 2},
  
  {460,  "CTEQ6.6A0",     "description", CTEQ66A_TBL_PATH, "ctq66.a0.pds",  2},
  {461,  "CTEQ6.6A1",     "description", CTEQ66A_TBL_PATH, "ctq66.a1.pds",  2},
  {462,  "CTEQ6.6A2",     "description", CTEQ66A_TBL_PATH, "ctq66.a2.pds",  2},
  {463,  "CTEQ6.6A3",     "description", CTEQ66A_TBL_PATH, "ctq66.a3.pds",  2},
  
  {450,  "CTEQ6.6C0",     "description", CTEQ66C_TBL_PATH, "ctq66.c0.pds",  2},
  {451,  "CTEQ6.6C1",     "description", CTEQ66C_TBL_PATH, "ctq66.c1.pds",  2},
  {452,  "CTEQ6.6C2",     "description", CTEQ66C_TBL_PATH, "ctq66.c2.pds",  2},
  {453,  "CTEQ6.6C3",     "description", CTEQ66C_TBL_PATH, "ctq66.c3.pds",  2},
  
  
  /*  End of the list */
  {0,0,0,0,0}
};



const cteq_pdfset_t *cteq_pdfset_database = __cteq_pdfset_database;


const cteq_pdfset_t * 
cteq_pdfset_find(const cteq_pdfset_t *pdflist, const char *name)
{
  while(pdflist->name) {
	if(strcmp(name, pdflist->name) == 0) return pdflist;
	++pdflist;
  }

  return 0;
}

const cteq_pdfset_t * 
cteq_pdfset_find_id(const cteq_pdfset_t *pdflist, int id)
{
  while(pdflist->name) {
	if(id == pdflist->id) return pdflist;
	++pdflist;
  }
  
  return 0;
}







