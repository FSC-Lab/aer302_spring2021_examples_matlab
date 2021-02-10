% plotprn.m - my personal plot print function
% 
% created on: 27-Feb-00
% updated on: 08-Sep-00
%
% Revision History
%
% 08-Sep-00: add saving as Bitmap file option

while (1)
   disp(' Print options: 1) no print 2) to printer 3) to PS file. 4) to BITMAP file');
   PrOp = input(' ??? (1,2,3,4) ');
   if (PrOp == 1)
      close;
      break;    
      elseif (PrOp == 2)
         %print  -dps;
         print
         close;
         break;
      elseif (PrOp == 3)
         PSfileName = input('input PS filename without extension: ','s');
         eval(['print ' PSfileName]);
         close;
         break;
      elseif (PrOp == 4)
         BITfileName = input('input BITMAP filename without extension: ','s');
         eval(['print -dbitmap ' BITfileName]);
         close;
         break;
      else
         disp(' Wrong answer! ');
      end % of if
   end % of while
   
