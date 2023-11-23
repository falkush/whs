package whs;

import java.awt.AWTException;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Cursor;
import java.awt.MouseInfo;
import java.awt.Point;
import java.awt.Robot;
import java.awt.Toolkit;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferByte;
import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.WindowConstants;

public class whs extends JPanel implements KeyListener, FocusListener {
	private static final long serialVersionUID = 1L;

	 static public char c;
	 static public boolean focus;
	 
	 static BufferedImage invcursimg = new BufferedImage(16, 16, BufferedImage.TYPE_INT_ARGB);
	 static Cursor invcurs = Toolkit.getDefaultToolkit().createCustomCursor(invcursimg, new Point(0, 0), "");
	 
	 static final JFrame frame = new JFrame("Press ESC to stop, click to resume");
	 static final JLabel label = new JLabel();
	 
	 static int centralx, centraly;
	 
	 static boolean holdw=false;
	 static boolean holda=false;
	 static boolean holds=false;
	 static boolean holdd=false;
	
	 static boolean holdz=false;
	 static boolean holdx=false;
	 static boolean holdc=false;
	 static boolean holdv=false;
	 static boolean resettot=false;
	 
	 static int n=3;
	

	 public static void main(String[] args) throws AWTException, InterruptedException {
		
		final int height=180;
		final int width=320;
		
		double dist=(double)1;
		double sqsz=(double)0.01;
		
		double bhsize=2*Math.PI;
		double whsize=Math.PI/2d;
		double roomsize=10d;
		double schecker=1d;
		
		
		int i,j,pw,ph;
		
		int[][] others = new int[3][2];
		
		double tmpl;
		double dc;
		double uhpy=1d;
		boolean flip=false;
		boolean flipg=false;
		double sx,sy,sz;
		double dotp,leangle;
		double[] vecperp = new double[3];
		double vecperpn;
		double lecos;
		double lesign;
		double lecossq;
		double lerayon;
		double lexw;
		double lex1,lex2;
		double dist2;
		double angx1;
		double distrem;
		double[] vectmpx = new double[3];
		double exitlgt,exitangle;
		double exitangle2;
		double[] pos2 = new double[3];
		double[] vecn2 = new double[3];
		
		double sxn,syn,szn;
		double theta,dist1,ang,ley2,newang;
		
		double dotp1,dotp2,projn,ang1;
			double[] newv = new double[3];
			double[] proj = new double[3];
		
		double tcont;
		double qa,qb,qc;
		double discr;
		double tmin,tmincoord;
		double tsol;
		double alpha;
		int setcoord;
		double[] coll = new double[3];
		int checker;
		int[] ctmp = new int[3];

		int currentpix;
		
		focus=true;
		
		final Robot robot = new Robot();
	    
	    
	    BufferedImage image = new BufferedImage(width,height,BufferedImage.TYPE_4BYTE_ABGR);

	    
		double msqsz=-sqsz;
		double multy=((double)(1-width))*sqsz/(double)2;
		double multz=((double)(height-1))*sqsz/(double)2;
		
		double[] addy = new double[n];
		double[] addz = new double[n];
		double[] vectmp = new double[n];
		double[] vecn = new double[n];
		double[][] vecl = new double[height][width];
		

			
		double[] pos = new double[n];
		
		
		int mousx=0;
		int mousy=0;
		
	
		
		double[][] x = new double[n][n];
		
		
		double[][] newx = new double[n][n];
		
		double[] vec = new double[n];

		double rsin = (double)Math.sin(Math.PI/32.0);
		double rcos = (double)Math.cos(Math.PI/32.0);
	
		others[0][0]=1;
		others[0][1]=2;
		others[1][0]=0;
		others[1][1]=2;
		others[2][0]=0;
		others[2][1]=1;

		for(i=0;i<n;i++) x[i][i]=1d;

		

		for(i=0;i<n;i++)
		{
			vec[i]=dist*x[0][i]+multy*x[1][i]+multz*x[2][i];
			addy[i]=sqsz*x[1][i];
			addz[i]=msqsz*x[2][i];
			pos[i]=0d;
		}
		pos[0]=-10d;
		
		
		for(ph=0;ph<height;ph++)
		{
			for(i=0;i<n;i++) vectmp[i]=vec[i];
			for(pw=0;pw<width;pw++)
			{
				vecl[ph][pw]=1d/veclgt(vec);
				for(i=0;i<n;i++) vec[i]+=addy[i];
			}
			for(i=0;i<n;i++) vec[i]=vectmp[i]+addz[i];
		}
	  
	    frame.setResizable(false);
	    frame.getContentPane().add(new whs());
	    frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
	    frame.setSize(width, height);
	    
	    frame.getContentPane().addMouseListener(new MouseAdapter() {            
	    	   @Override
	    	   public void mouseClicked(MouseEvent e) {
	    		  
	    	      focus=true;
	    	      centralx=frame.getX()+(width/2);
	    	      centraly=frame.getY()+(height/2);
	    	      frame.getContentPane().setCursor(invcurs);
	    	   }
	    	});
	    
	    frame.getContentPane().setCursor(invcurs);
	    frame.getContentPane().setBackground( Color.BLACK );
	    
	    final byte[] pixels =((DataBufferByte) image.getRaster().getDataBuffer()).getData();

       label.setIcon(new ImageIcon(image));
       frame.getContentPane().add(label,BorderLayout.CENTER);
       frame.setLocationRelativeTo(null);
       frame.pack();
       frame.setVisible(true);
	   
       centralx=frame.getX()+(width/2);
       centraly=frame.getY()+(height/2);
       
       
	    while(true)
	    {
	    	Thread.sleep(10);
	    	
	    	if(focus)
	    	{
	    			if(holdz) bhsize+=0.1d;
	    			if(holdx && bhsize>0.1d) bhsize-=0.1d;
	    			if(holdc) whsize+=0.1d;
	    			if(holdv && whsize>0.1d) whsize-=0.1d;
	    			
	    			
			    	
			    	if(resettot)
			    	{
			    		for(i=0;i<n;i++) {pos[i]=0d; for(j=0;j<n;j++) if(i==j) x[i][j]=1; else x[i][j]=0;}
			    		pos[0]=-5d;
			    		bhsize=2*Math.PI;
			    		whsize=Math.PI/2d;
			    		resettot=false;
			    	}
		    		if(holdw)
		    		{
		    			
		    			dc=Math.sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
				    	uhpy=1d/dc;
				    	
		    			if(dc>1d) {
		    			 for(i=0;i<n;i++) pos[i]+=x[0][i]/10.0d;
		    			}
		    			else
		    			{
		    				
		    				sxn=-pos[0]/dc;
		    				syn=-pos[1]/dc;
		    				szn=-pos[2]/dc;
		    				
		    				dotp=sxn*x[0][0]+syn*x[0][1]+szn*x[0][2];
		    				leangle=Math.acos(dotp);
		    				
		    				theta=Math.PI/2d-leangle;
		    					
			    				dist1=arsinh(Math.tan(theta));
			    				ang=Math.atan(Math.sinh(0.1d-dist1));
			    				
			    				lerayon=uhpy/Math.sin(Math.PI/2d+theta);
			    				
			    				lex1=lerayon*Math.cos(Math.PI/2d+theta);
			    				
			    				lex2=lerayon*Math.cos(Math.PI/2d-ang);
			    				ley2=lerayon*Math.sin(Math.PI/2d-ang);
			    				
			    				sxn=-sxn;
			    				syn=-syn;
			    				szn=-szn;
			    				dotp=-dotp;
			    				leangle=Math.acos(dotp);
			    				
			    				vecperp[0]=x[0][0]-dotp*sxn;
			 	    			vecperp[1]=x[0][1]-dotp*syn;
			 	    			vecperp[2]=x[0][2]-dotp*szn;
			 	    				
			 	    			vecperpn = Math.sqrt(vecperp[0]*vecperp[0]+vecperp[1]*vecperp[1]+vecperp[2]*vecperp[2]);
			 	    				
			 	    			vecperp[0]/=vecperpn;
			 	    			vecperp[1]/=vecperpn;
			 	    			vecperp[2]/=vecperpn;
			 	    			
			 	    			
			 	    			
			 	    			if(ley2<whsize)
			 	    			{
			 	    			
				 	    				exitlgt=lex2-lex1;
					 	    			exitangle=exitlgt%bhsize;
					 	    			exitangle=(2d*Math.PI/bhsize)*exitangle;
					 	    			newang=Math.PI/2d+ang;
					 	    			
					 	    			pos[0]=Math.cos(exitangle)*sxn+Math.sin(exitangle)*vecperp[0];
					 	    			pos[1]=Math.cos(exitangle)*syn+Math.sin(exitangle)*vecperp[1];
					 	    			pos[2]=Math.cos(exitangle)*szn+Math.sin(exitangle)*vecperp[2];
					 	    			
					 	    			pos[0]*=1/ley2;
					 	    			pos[1]*=1/ley2;
					 	    			pos[2]*=1/ley2;	
				 	    				
				 	    			newv[0]=Math.cos(Math.PI+exitangle-newang)*sxn+Math.sin(Math.PI+exitangle-newang)*vecperp[0];
				 	    			newv[1]=Math.cos(Math.PI+exitangle-newang)*syn+Math.sin(Math.PI+exitangle-newang)*vecperp[1];
				 	    			newv[2]=Math.cos(Math.PI+exitangle-newang)*szn+Math.sin(Math.PI+exitangle-newang)*vecperp[2];
				 	    			
				 	    			
				 	    			for(j=1;j<3;j++) 
				 	    			{
					 	    			dotp1=newv[0]*x[0][0]+newv[1]*x[0][1]+newv[2]*x[0][2];
					 	    			dotp2=newv[0]*x[j][0]+newv[1]*x[j][1]+newv[2]*x[j][2];
					 	    			
					 	    			proj[0]=x[0][0]*dotp1+x[j][0]*dotp2;
					 	    			proj[1]=x[0][1]*dotp1+x[j][1]*dotp2;
					 	    			proj[2]=x[0][2]*dotp1+x[j][2]*dotp2;
					 	    			
					 	    			projn = Math.sqrt(proj[0]*proj[0]+proj[1]*proj[1]+proj[2]*proj[2]);
				 	    				
					 	    			proj[0]/=projn;
					 	    			proj[1]/=projn;
					 	    			proj[2]/=projn;
					 	    			
					 	    			tmpl=proj[0]*x[0][0]+proj[1]*x[0][1]+proj[2]*x[0][2];
					 	    			if(Math.abs(tmpl)>1d) ang1=0;
					 	    			else
					 	    			{
						 	    			ang1=Math.acos(tmpl);
						 	    			ang1*=Math.signum(proj[0]*x[j][0]+proj[1]*x[j][1]+proj[2]*x[j][2]);
					 	    			}
					 	    			for(i=0;i<3;i++)
							    		{
							    			newx[0][i]=x[0][i]*Math.cos(ang1)+Math.sin(ang1)*x[j][i];
							    			newx[j][i]=x[j][i]*Math.cos(ang1)-Math.sin(ang1)*x[0][i];
							    			x[0][i]=newx[0][i];
							    			x[j][i]=newx[j][i];
							    		}
					 	    			
				 	    			}
			 	    			}
			 	    			else
			 	    			{
			 	    				
			 	    				flipg=!flipg;
	
			 	    				lex2=-Math.sqrt(lerayon*lerayon-whsize*whsize);
			 	    				
			 	    				angx1=Math.PI/2d - Math.asin(whsize/lerayon);
				 	    			dist2=arsinh(Math.tan(angx1));
				 	    			distrem=0.1d-dist1+dist2;
				 	    			ang=Math.atan(Math.sinh(0.1d-distrem-dist1));
				 	    			
				 	    			exitlgt=lex2-lex1;
				 	    			exitangle=exitlgt%bhsize;
				 	    			exitangle=(2d*Math.PI/bhsize)*exitangle;
				 	    			newang=Math.PI/2d+ang;
				 	    			
				 	    			pos[0]=Math.cos(exitangle)*sxn+Math.sin(exitangle)*vecperp[0];
				 	    			pos[1]=Math.cos(exitangle)*syn+Math.sin(exitangle)*vecperp[1];
				 	    			pos[2]=Math.cos(exitangle)*szn+Math.sin(exitangle)*vecperp[2];
				 	    			
				 	    			pos[0]*=1/whsize;
				 	    			pos[1]*=1/whsize;
				 	    			pos[2]*=1/whsize;	
			 	    				
			 	    			newv[0]=Math.cos(Math.PI+exitangle-newang)*sxn+Math.sin(Math.PI+exitangle-newang)*vecperp[0];
			 	    			newv[1]=Math.cos(Math.PI+exitangle-newang)*syn+Math.sin(Math.PI+exitangle-newang)*vecperp[1];
			 	    			newv[2]=Math.cos(Math.PI+exitangle-newang)*szn+Math.sin(Math.PI+exitangle-newang)*vecperp[2];
			 	    			
			 	    			
			 	    			for(j=1;j<3;j++) 
			 	    			{
				 	    			dotp1=newv[0]*x[0][0]+newv[1]*x[0][1]+newv[2]*x[0][2];
				 	    			dotp2=newv[0]*x[j][0]+newv[1]*x[j][1]+newv[2]*x[j][2];
				 	    			
				 	    			proj[0]=x[0][0]*dotp1+x[j][0]*dotp2;
				 	    			proj[1]=x[0][1]*dotp1+x[j][1]*dotp2;
				 	    			proj[2]=x[0][2]*dotp1+x[j][2]*dotp2;
				 	    			
				 	    			projn = Math.sqrt(proj[0]*proj[0]+proj[1]*proj[1]+proj[2]*proj[2]);
			 	    				
				 	    			proj[0]/=projn;
				 	    			proj[1]/=projn;
				 	    			proj[2]/=projn;
				 	    			
				 	    			tmpl=proj[0]*x[0][0]+proj[1]*x[0][1]+proj[2]*x[0][2];
				 	    			if(Math.abs(tmpl)>1d) ang1=0;
				 	    			else
				 	    			{
					 	    			ang1=Math.acos(tmpl);
					 	    			ang1*=Math.signum(proj[0]*x[j][0]+proj[1]*x[j][1]+proj[2]*x[j][2]);
				 	    			}
				 	    			for(i=0;i<3;i++)
						    		{
						    			newx[0][i]=x[0][i]*Math.cos(ang1)+Math.sin(ang1)*x[j][i];
						    			newx[j][i]=x[j][i]*Math.cos(ang1)-Math.sin(ang1)*x[0][i];
						    			x[0][i]=newx[0][i];
						    			x[j][i]=newx[j][i];
						    		}
				 	    			
			 	    			}
				 	    			
				 	    			
				 	    			
				 	    			
			 	    			dc=Math.sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
						    	uhpy=1d/dc;
						    	sxn=pos[0]/dc;
			    				syn=pos[1]/dc;
			    				szn=pos[2]/dc;
			    				
			    				for(j=0;j<3;j++) {
				    				dotp=sxn*x[j][0]+syn*x[j][1]+szn*x[j][2];	
				    				
				    				vectmpx[0]=sxn*dotp;
				    				vectmpx[1]=syn*dotp;
				    				vectmpx[2]=szn*dotp;
				    				
				    				x[j][0]-=2d*vectmpx[0];
				    				x[j][1]-=2d*vectmpx[1];
				    				x[j][2]-=2d*vectmpx[2];
			    				}
			    				
			    				sxn=-pos[0]/dc;
			    				syn=-pos[1]/dc;
			    				szn=-pos[2]/dc;
			    				
			    				dotp=sxn*x[0][0]+syn*x[0][1]+szn*x[0][2];
			    				leangle=Math.acos(dotp);
			    				
			    				theta=Math.PI/2d-leangle;
			    					
				    				dist1=arsinh(Math.tan(theta));
				    				ang=Math.atan(Math.sinh(distrem-dist1));
				    				
				    				lerayon=uhpy/Math.sin(Math.PI/2d+theta);
				    				
				    				lex1=lerayon*Math.cos(Math.PI/2d+theta);
				    				
				    				lex2=lerayon*Math.cos(Math.PI/2d-ang);
				    				ley2=lerayon*Math.sin(Math.PI/2d-ang);
				    				
				    				sxn=-sxn;
				    				syn=-syn;
				    				szn=-szn;
				    				dotp=-dotp;
				    				leangle=Math.acos(dotp);
				    				
				    				vecperp[0]=x[0][0]-dotp*sxn;
				 	    			vecperp[1]=x[0][1]-dotp*syn;
				 	    			vecperp[2]=x[0][2]-dotp*szn;
				 	    				
				 	    			vecperpn = Math.sqrt(vecperp[0]*vecperp[0]+vecperp[1]*vecperp[1]+vecperp[2]*vecperp[2]);
				 	    				
				 	    			vecperp[0]/=vecperpn;
				 	    			vecperp[1]/=vecperpn;
				 	    			vecperp[2]/=vecperpn;

				 	    			
					 	    				exitlgt=lex2-lex1;
						 	    			exitangle=exitlgt%bhsize;
						 	    			exitangle=(2d*Math.PI/bhsize)*exitangle;
						 	    			newang=Math.PI/2d+ang;
						 	    			
						 	    			pos[0]=Math.cos(exitangle)*sxn+Math.sin(exitangle)*vecperp[0];
						 	    			pos[1]=Math.cos(exitangle)*syn+Math.sin(exitangle)*vecperp[1];
						 	    			pos[2]=Math.cos(exitangle)*szn+Math.sin(exitangle)*vecperp[2];
						 	    			
						 	    			pos[0]*=1/ley2;
						 	    			pos[1]*=1/ley2;
						 	    			pos[2]*=1/ley2;	
					 	    				
					 	    			newv[0]=Math.cos(Math.PI+exitangle-newang)*sxn+Math.sin(Math.PI+exitangle-newang)*vecperp[0];
					 	    			newv[1]=Math.cos(Math.PI+exitangle-newang)*syn+Math.sin(Math.PI+exitangle-newang)*vecperp[1];
					 	    			newv[2]=Math.cos(Math.PI+exitangle-newang)*szn+Math.sin(Math.PI+exitangle-newang)*vecperp[2];
					 	    			
					 	    			
					 	    			for(j=1;j<3;j++) 
					 	    			{
						 	    			dotp1=newv[0]*x[0][0]+newv[1]*x[0][1]+newv[2]*x[0][2];
						 	    			dotp2=newv[0]*x[j][0]+newv[1]*x[j][1]+newv[2]*x[j][2];
						 	    			
						 	    			proj[0]=x[0][0]*dotp1+x[j][0]*dotp2;
						 	    			proj[1]=x[0][1]*dotp1+x[j][1]*dotp2;
						 	    			proj[2]=x[0][2]*dotp1+x[j][2]*dotp2;
						 	    			
						 	    			projn = Math.sqrt(proj[0]*proj[0]+proj[1]*proj[1]+proj[2]*proj[2]);
					 	    				
						 	    			proj[0]/=projn;
						 	    			proj[1]/=projn;
						 	    			proj[2]/=projn;
						 	    			
						 	    			tmpl=proj[0]*x[0][0]+proj[1]*x[0][1]+proj[2]*x[0][2];
						 	    			if(Math.abs(tmpl)>1d) ang1=0;
						 	    			else
						 	    			{
							 	    			ang1=Math.acos(tmpl);
							 	    			ang1*=Math.signum(proj[0]*x[j][0]+proj[1]*x[j][1]+proj[2]*x[j][2]);
						 	    			}
						 	    			for(i=0;i<3;i++)
								    		{
								    			newx[0][i]=x[0][i]*Math.cos(ang1)+Math.sin(ang1)*x[j][i];
								    			newx[j][i]=x[j][i]*Math.cos(ang1)-Math.sin(ang1)*x[0][i];
								    			x[0][i]=newx[0][i];
								    			x[j][i]=newx[j][i];
								    		}
						 	    			
					 	    			}
			 	    			}
			 	    			
			 	    			
		    				}
		    			
		    			 
		    			 
		    			 
		    		}
		    		if(holds)
		    		{
		    			dc=Math.sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
				    	uhpy=1d/dc;
				    	
		    			if(dc>1d) {
		    			 for(i=0;i<n;i++) pos[i]-=x[0][i]/10.0d;
		    			}
		    			else
		    			{
		    				x[0][0]*=-1d;
		    				x[0][1]*=-1d;
		    				x[0][2]*=-1d;
		    				
		    				sxn=-pos[0]/dc;
		    				syn=-pos[1]/dc;
		    				szn=-pos[2]/dc;
		    				
		    				dotp=sxn*x[0][0]+syn*x[0][1]+szn*x[0][2];
		    				leangle=Math.acos(dotp);
		    				
		    				theta=Math.PI/2d-leangle;
		    					
			    				dist1=arsinh(Math.tan(theta));
			    				ang=Math.atan(Math.sinh(0.1d-dist1));
			    				
			    				lerayon=uhpy/Math.sin(Math.PI/2d+theta);
			    				
			    				lex1=lerayon*Math.cos(Math.PI/2d+theta);
			    				
			    				lex2=lerayon*Math.cos(Math.PI/2d-ang);
			    				ley2=lerayon*Math.sin(Math.PI/2d-ang);
			    				
			    				sxn=-sxn;
			    				syn=-syn;
			    				szn=-szn;
			    				dotp=-dotp;
			    				leangle=Math.acos(dotp);
			    				
			    				vecperp[0]=x[0][0]-dotp*sxn;
			 	    			vecperp[1]=x[0][1]-dotp*syn;
			 	    			vecperp[2]=x[0][2]-dotp*szn;
			 	    				
			 	    			vecperpn = Math.sqrt(vecperp[0]*vecperp[0]+vecperp[1]*vecperp[1]+vecperp[2]*vecperp[2]);
			 	    				
			 	    			vecperp[0]/=vecperpn;
			 	    			vecperp[1]/=vecperpn;
			 	    			vecperp[2]/=vecperpn;
			 	    			
			 	    			
			 	    			
			 	    			if(ley2<whsize)
			 	    			{
			 	    			
				 	    				exitlgt=lex2-lex1;
					 	    			exitangle=exitlgt%bhsize;
					 	    			exitangle=(2d*Math.PI/bhsize)*exitangle;
					 	    			newang=Math.PI/2d+ang;
					 	    			
					 	    			pos[0]=Math.cos(exitangle)*sxn+Math.sin(exitangle)*vecperp[0];
					 	    			pos[1]=Math.cos(exitangle)*syn+Math.sin(exitangle)*vecperp[1];
					 	    			pos[2]=Math.cos(exitangle)*szn+Math.sin(exitangle)*vecperp[2];
					 	    			
					 	    			pos[0]*=1/ley2;
					 	    			pos[1]*=1/ley2;
					 	    			pos[2]*=1/ley2;	
				 	    				
				 	    			newv[0]=Math.cos(Math.PI+exitangle-newang)*sxn+Math.sin(Math.PI+exitangle-newang)*vecperp[0];
				 	    			newv[1]=Math.cos(Math.PI+exitangle-newang)*syn+Math.sin(Math.PI+exitangle-newang)*vecperp[1];
				 	    			newv[2]=Math.cos(Math.PI+exitangle-newang)*szn+Math.sin(Math.PI+exitangle-newang)*vecperp[2];
				 	    			
				 	    			
				 	    			for(j=1;j<3;j++) 
				 	    			{
					 	    			dotp1=newv[0]*x[0][0]+newv[1]*x[0][1]+newv[2]*x[0][2];
					 	    			dotp2=newv[0]*x[j][0]+newv[1]*x[j][1]+newv[2]*x[j][2];
					 	    			
					 	    			proj[0]=x[0][0]*dotp1+x[j][0]*dotp2;
					 	    			proj[1]=x[0][1]*dotp1+x[j][1]*dotp2;
					 	    			proj[2]=x[0][2]*dotp1+x[j][2]*dotp2;
					 	    			
					 	    			projn = Math.sqrt(proj[0]*proj[0]+proj[1]*proj[1]+proj[2]*proj[2]);
				 	    				
					 	    			proj[0]/=projn;
					 	    			proj[1]/=projn;
					 	    			proj[2]/=projn;
					 	    			
					 	    			tmpl=proj[0]*x[0][0]+proj[1]*x[0][1]+proj[2]*x[0][2];
					 	    			if(Math.abs(tmpl)>1d) ang1=0;
					 	    			else
					 	    			{
						 	    			ang1=Math.acos(tmpl);
						 	    			ang1*=Math.signum(proj[0]*x[j][0]+proj[1]*x[j][1]+proj[2]*x[j][2]);
					 	    			}
					 	    			for(i=0;i<3;i++)
							    		{
							    			newx[0][i]=x[0][i]*Math.cos(ang1)+Math.sin(ang1)*x[j][i];
							    			newx[j][i]=x[j][i]*Math.cos(ang1)-Math.sin(ang1)*x[0][i];
							    			x[0][i]=newx[0][i];
							    			x[j][i]=newx[j][i];
							    		}
					 	    			
				 	    			}
			 	    			}
			 	    			else
			 	    			{
			 	    				
			 	    				flipg=!flipg;
	
			 	    				lex2=-Math.sqrt(lerayon*lerayon-whsize*whsize);
			 	    				
			 	    				angx1=Math.PI/2d - Math.asin(whsize/lerayon);
				 	    			dist2=arsinh(Math.tan(angx1));
				 	    			distrem=0.1d-dist1+dist2;
				 	    			ang=Math.atan(Math.sinh(0.1d-distrem-dist1));
				 	    			
				 	    			exitlgt=lex2-lex1;
				 	    			exitangle=exitlgt%bhsize;
				 	    			exitangle=(2d*Math.PI/bhsize)*exitangle;
				 	    			newang=Math.PI/2d+ang;
				 	    			
				 	    			pos[0]=Math.cos(exitangle)*sxn+Math.sin(exitangle)*vecperp[0];
				 	    			pos[1]=Math.cos(exitangle)*syn+Math.sin(exitangle)*vecperp[1];
				 	    			pos[2]=Math.cos(exitangle)*szn+Math.sin(exitangle)*vecperp[2];
				 	    			
				 	    			pos[0]*=1/whsize;
				 	    			pos[1]*=1/whsize;
				 	    			pos[2]*=1/whsize;	
			 	    				
			 	    			newv[0]=Math.cos(Math.PI+exitangle-newang)*sxn+Math.sin(Math.PI+exitangle-newang)*vecperp[0];
			 	    			newv[1]=Math.cos(Math.PI+exitangle-newang)*syn+Math.sin(Math.PI+exitangle-newang)*vecperp[1];
			 	    			newv[2]=Math.cos(Math.PI+exitangle-newang)*szn+Math.sin(Math.PI+exitangle-newang)*vecperp[2];
			 	    			
			 	    			
			 	    			for(j=1;j<3;j++) 
			 	    			{
				 	    			dotp1=newv[0]*x[0][0]+newv[1]*x[0][1]+newv[2]*x[0][2];
				 	    			dotp2=newv[0]*x[j][0]+newv[1]*x[j][1]+newv[2]*x[j][2];
				 	    			
				 	    			proj[0]=x[0][0]*dotp1+x[j][0]*dotp2;
				 	    			proj[1]=x[0][1]*dotp1+x[j][1]*dotp2;
				 	    			proj[2]=x[0][2]*dotp1+x[j][2]*dotp2;
				 	    			
				 	    			projn = Math.sqrt(proj[0]*proj[0]+proj[1]*proj[1]+proj[2]*proj[2]);
			 	    				
				 	    			proj[0]/=projn;
				 	    			proj[1]/=projn;
				 	    			proj[2]/=projn;
				 	    			
				 	    			tmpl=proj[0]*x[0][0]+proj[1]*x[0][1]+proj[2]*x[0][2];
				 	    			if(Math.abs(tmpl)>1d) ang1=0;
				 	    			else
				 	    			{
					 	    			ang1=Math.acos(tmpl);
					 	    			ang1*=Math.signum(proj[0]*x[j][0]+proj[1]*x[j][1]+proj[2]*x[j][2]);
				 	    			}
				 	    			for(i=0;i<3;i++)
						    		{
						    			newx[0][i]=x[0][i]*Math.cos(ang1)+Math.sin(ang1)*x[j][i];
						    			newx[j][i]=x[j][i]*Math.cos(ang1)-Math.sin(ang1)*x[0][i];
						    			x[0][i]=newx[0][i];
						    			x[j][i]=newx[j][i];
						    		}
				 	    			
			 	    			}
				 	    			
				 	    			
				 	    			
				 	    			
			 	    			dc=Math.sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
						    	uhpy=1d/dc;
						    	sxn=pos[0]/dc;
			    				syn=pos[1]/dc;
			    				szn=pos[2]/dc;
			    				
			    				for(j=0;j<3;j++) {
				    				dotp=sxn*x[j][0]+syn*x[j][1]+szn*x[j][2];	
				    				
				    				vectmpx[0]=sxn*dotp;
				    				vectmpx[1]=syn*dotp;
				    				vectmpx[2]=szn*dotp;
				    				
				    				x[j][0]-=2d*vectmpx[0];
				    				x[j][1]-=2d*vectmpx[1];
				    				x[j][2]-=2d*vectmpx[2];
			    				}
			    				
			    				sxn=-pos[0]/dc;
			    				syn=-pos[1]/dc;
			    				szn=-pos[2]/dc;
			    				
			    				dotp=sxn*x[0][0]+syn*x[0][1]+szn*x[0][2];
			    				leangle=Math.acos(dotp);
			    				
			    				theta=Math.PI/2d-leangle;
			    					
				    				dist1=arsinh(Math.tan(theta));
				    				ang=Math.atan(Math.sinh(distrem-dist1));
				    				
				    				lerayon=uhpy/Math.sin(Math.PI/2d+theta);
				    				
				    				lex1=lerayon*Math.cos(Math.PI/2d+theta);
				    				
				    				lex2=lerayon*Math.cos(Math.PI/2d-ang);
				    				ley2=lerayon*Math.sin(Math.PI/2d-ang);
				    				
				    				sxn=-sxn;
				    				syn=-syn;
				    				szn=-szn;
				    				dotp=-dotp;
				    				leangle=Math.acos(dotp);
				    				
				    				vecperp[0]=x[0][0]-dotp*sxn;
				 	    			vecperp[1]=x[0][1]-dotp*syn;
				 	    			vecperp[2]=x[0][2]-dotp*szn;
				 	    				
				 	    			vecperpn = Math.sqrt(vecperp[0]*vecperp[0]+vecperp[1]*vecperp[1]+vecperp[2]*vecperp[2]);
				 	    				
				 	    			vecperp[0]/=vecperpn;
				 	    			vecperp[1]/=vecperpn;
				 	    			vecperp[2]/=vecperpn;

				 	    			
					 	    				exitlgt=lex2-lex1;
						 	    			exitangle=exitlgt%bhsize;
						 	    			exitangle=(2d*Math.PI/bhsize)*exitangle;
						 	    			newang=Math.PI/2d+ang;
						 	    			
						 	    			pos[0]=Math.cos(exitangle)*sxn+Math.sin(exitangle)*vecperp[0];
						 	    			pos[1]=Math.cos(exitangle)*syn+Math.sin(exitangle)*vecperp[1];
						 	    			pos[2]=Math.cos(exitangle)*szn+Math.sin(exitangle)*vecperp[2];
						 	    			
						 	    			pos[0]*=1/ley2;
						 	    			pos[1]*=1/ley2;
						 	    			pos[2]*=1/ley2;	
					 	    				
					 	    			newv[0]=Math.cos(Math.PI+exitangle-newang)*sxn+Math.sin(Math.PI+exitangle-newang)*vecperp[0];
					 	    			newv[1]=Math.cos(Math.PI+exitangle-newang)*syn+Math.sin(Math.PI+exitangle-newang)*vecperp[1];
					 	    			newv[2]=Math.cos(Math.PI+exitangle-newang)*szn+Math.sin(Math.PI+exitangle-newang)*vecperp[2];
					 	    			
					 	    			
					 	    			for(j=1;j<3;j++) 
					 	    			{
						 	    			dotp1=newv[0]*x[0][0]+newv[1]*x[0][1]+newv[2]*x[0][2];
						 	    			dotp2=newv[0]*x[j][0]+newv[1]*x[j][1]+newv[2]*x[j][2];
						 	    			
						 	    			proj[0]=x[0][0]*dotp1+x[j][0]*dotp2;
						 	    			proj[1]=x[0][1]*dotp1+x[j][1]*dotp2;
						 	    			proj[2]=x[0][2]*dotp1+x[j][2]*dotp2;
						 	    			
						 	    			projn = Math.sqrt(proj[0]*proj[0]+proj[1]*proj[1]+proj[2]*proj[2]);
					 	    				
						 	    			proj[0]/=projn;
						 	    			proj[1]/=projn;
						 	    			proj[2]/=projn;
						 	    			
						 	    			tmpl=proj[0]*x[0][0]+proj[1]*x[0][1]+proj[2]*x[0][2];
						 	    			if(Math.abs(tmpl)>1d) ang1=0;
						 	    			else
						 	    			{
							 	    			ang1=Math.acos(tmpl);
							 	    			ang1*=Math.signum(proj[0]*x[j][0]+proj[1]*x[j][1]+proj[2]*x[j][2]);
						 	    			}
						 	    			for(i=0;i<3;i++)
								    		{
								    			newx[0][i]=x[0][i]*Math.cos(ang1)+Math.sin(ang1)*x[j][i];
								    			newx[j][i]=x[j][i]*Math.cos(ang1)-Math.sin(ang1)*x[0][i];
								    			x[0][i]=newx[0][i];
								    			x[j][i]=newx[j][i];
								    		}
						 	    			
					 	    			}
			 	    			}
			 	    			
			 	    			x[0][0]*=-1d;
			    				x[0][1]*=-1d;
			    				x[0][2]*=-1d;
		    				}

		    		}
		    		if(holda)
		    		{
		    			dc=Math.sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
				    	uhpy=1d/dc;
				    	
		    			if(dc>1d) {
		    			 for(i=0;i<n;i++) pos[i]-=x[1][i]/10.0d;
		    			}
		    			else
		    			{
		    				
		    				x[1][0]*=-1d;
		    				x[1][1]*=-1d;
		    				x[1][2]*=-1d;
		    				
		    				sxn=-pos[0]/dc;
		    				syn=-pos[1]/dc;
		    				szn=-pos[2]/dc;
		    				
		    				dotp=sxn*x[1][0]+syn*x[1][1]+szn*x[1][2];
		    				leangle=Math.acos(dotp);
		    				
		    				theta=Math.PI/2d-leangle;
		    					
			    				dist1=arsinh(Math.tan(theta));
			    				ang=Math.atan(Math.sinh(0.1d-dist1));
			    				
			    				lerayon=uhpy/Math.sin(Math.PI/2d+theta);
			    				
			    				lex1=lerayon*Math.cos(Math.PI/2d+theta);
			    				
			    				lex2=lerayon*Math.cos(Math.PI/2d-ang);
			    				ley2=lerayon*Math.sin(Math.PI/2d-ang);
			    				
			    				sxn=-sxn;
			    				syn=-syn;
			    				szn=-szn;
			    				dotp=-dotp;
			    				leangle=Math.acos(dotp);
			    				
			    				vecperp[0]=x[1][0]-dotp*sxn;
			 	    			vecperp[1]=x[1][1]-dotp*syn;
			 	    			vecperp[2]=x[1][2]-dotp*szn;
			 	    				
			 	    			vecperpn = Math.sqrt(vecperp[0]*vecperp[0]+vecperp[1]*vecperp[1]+vecperp[2]*vecperp[2]);
			 	    				
			 	    			vecperp[0]/=vecperpn;
			 	    			vecperp[1]/=vecperpn;
			 	    			vecperp[2]/=vecperpn;
			 	    			
			 	    			
			 	    			
			 	    			if(ley2<whsize)
			 	    			{
			 	    			
				 	    				exitlgt=lex2-lex1;
					 	    			exitangle=exitlgt%bhsize;
					 	    			exitangle=(2d*Math.PI/bhsize)*exitangle;
					 	    			newang=Math.PI/2d+ang;
					 	    			
					 	    			pos[0]=Math.cos(exitangle)*sxn+Math.sin(exitangle)*vecperp[0];
					 	    			pos[1]=Math.cos(exitangle)*syn+Math.sin(exitangle)*vecperp[1];
					 	    			pos[2]=Math.cos(exitangle)*szn+Math.sin(exitangle)*vecperp[2];
					 	    			
					 	    			pos[0]*=1/ley2;
					 	    			pos[1]*=1/ley2;
					 	    			pos[2]*=1/ley2;	
				 	    				
				 	    			newv[0]=Math.cos(Math.PI+exitangle-newang)*sxn+Math.sin(Math.PI+exitangle-newang)*vecperp[0];
				 	    			newv[1]=Math.cos(Math.PI+exitangle-newang)*syn+Math.sin(Math.PI+exitangle-newang)*vecperp[1];
				 	    			newv[2]=Math.cos(Math.PI+exitangle-newang)*szn+Math.sin(Math.PI+exitangle-newang)*vecperp[2];
				 	    			
				 	    			
				 	    			for(j=0;j<2;j++) 
				 	    			{
					 	    			dotp1=newv[0]*x[1][0]+newv[1]*x[1][1]+newv[2]*x[1][2];
					 	    			dotp2=newv[0]*x[others[1][j]][0]+newv[1]*x[others[1][j]][1]+newv[2]*x[others[1][j]][2];
					 	    			
					 	    			proj[0]=x[1][0]*dotp1+x[others[1][j]][0]*dotp2;
					 	    			proj[1]=x[1][1]*dotp1+x[others[1][j]][1]*dotp2;
					 	    			proj[2]=x[1][2]*dotp1+x[others[1][j]][2]*dotp2;
					 	    			
					 	    			projn = Math.sqrt(proj[0]*proj[0]+proj[1]*proj[1]+proj[2]*proj[2]);
				 	    				
					 	    			proj[0]/=projn;
					 	    			proj[1]/=projn;
					 	    			proj[2]/=projn;
					 	    			
					 	    			tmpl=proj[0]*x[1][0]+proj[1]*x[1][1]+proj[2]*x[1][2];
					 	    			if(Math.abs(tmpl)>1d) ang1=0;
					 	    			else
					 	    			{
						 	    			ang1=Math.acos(tmpl);
						 	    			ang1*=Math.signum(proj[0]*x[others[1][j]][0]+proj[1]*x[others[1][j]][1]+proj[2]*x[others[1][j]][2]);
					 	    			}
					 	    			for(i=0;i<3;i++)
							    		{
							    			newx[1][i]=x[1][i]*Math.cos(ang1)+Math.sin(ang1)*x[others[1][j]][i];
							    			newx[others[1][j]][i]=x[others[1][j]][i]*Math.cos(ang1)-Math.sin(ang1)*x[1][i];
							    			x[1][i]=newx[1][i];
							    			x[others[1][j]][i]=newx[others[1][j]][i];
							    		}
					 	    			
				 	    			}
			 	    			}
			 	    			else
			 	    			{
			 	    				
			 	    				flipg=!flipg;
	
			 	    				lex2=-Math.sqrt(lerayon*lerayon-whsize*whsize);
			 	    				
			 	    				angx1=Math.PI/2d - Math.asin(whsize/lerayon);
				 	    			dist2=arsinh(Math.tan(angx1));
				 	    			distrem=0.1d-dist1+dist2;
				 	    			ang=Math.atan(Math.sinh(0.1d-distrem-dist1));
				 	    			
				 	    			exitlgt=lex2-lex1;
				 	    			exitangle=exitlgt%bhsize;
				 	    			exitangle=(2d*Math.PI/bhsize)*exitangle;
				 	    			newang=Math.PI/2d+ang;
				 	    			
				 	    			pos[0]=Math.cos(exitangle)*sxn+Math.sin(exitangle)*vecperp[0];
				 	    			pos[1]=Math.cos(exitangle)*syn+Math.sin(exitangle)*vecperp[1];
				 	    			pos[2]=Math.cos(exitangle)*szn+Math.sin(exitangle)*vecperp[2];
				 	    			
				 	    			pos[0]*=1/whsize;
				 	    			pos[1]*=1/whsize;
				 	    			pos[2]*=1/whsize;	
			 	    				
			 	    			newv[0]=Math.cos(Math.PI+exitangle-newang)*sxn+Math.sin(Math.PI+exitangle-newang)*vecperp[0];
			 	    			newv[1]=Math.cos(Math.PI+exitangle-newang)*syn+Math.sin(Math.PI+exitangle-newang)*vecperp[1];
			 	    			newv[2]=Math.cos(Math.PI+exitangle-newang)*szn+Math.sin(Math.PI+exitangle-newang)*vecperp[2];
			 	    			
			 	    			
			 	    			for(j=0;j<2;j++) 
			 	    			{
				 	    			dotp1=newv[0]*x[1][0]+newv[1]*x[1][1]+newv[2]*x[1][2];
				 	    			dotp2=newv[0]*x[others[1][j]][0]+newv[1]*x[others[1][j]][1]+newv[2]*x[others[1][j]][2];
				 	    			
				 	    			proj[0]=x[1][0]*dotp1+x[others[1][j]][0]*dotp2;
				 	    			proj[1]=x[1][1]*dotp1+x[others[1][j]][1]*dotp2;
				 	    			proj[2]=x[1][2]*dotp1+x[others[1][j]][2]*dotp2;
				 	    			
				 	    			projn = Math.sqrt(proj[0]*proj[0]+proj[1]*proj[1]+proj[2]*proj[2]);
			 	    				
				 	    			proj[0]/=projn;
				 	    			proj[1]/=projn;
				 	    			proj[2]/=projn;
				 	    			
				 	    			tmpl=proj[0]*x[1][0]+proj[1]*x[1][1]+proj[2]*x[1][2];
				 	    			if(Math.abs(tmpl)>1d) ang1=0;
				 	    			else
				 	    			{
					 	    			ang1=Math.acos(tmpl);
					 	    			ang1*=Math.signum(proj[0]*x[others[1][j]][0]+proj[1]*x[others[1][j]][1]+proj[2]*x[others[1][j]][2]);
				 	    			}
				 	    			for(i=0;i<3;i++)
						    		{
						    			newx[1][i]=x[1][i]*Math.cos(ang1)+Math.sin(ang1)*x[others[1][j]][i];
						    			newx[others[1][j]][i]=x[others[1][j]][i]*Math.cos(ang1)-Math.sin(ang1)*x[1][i];
						    			x[1][i]=newx[1][i];
						    			x[others[1][j]][i]=newx[others[1][j]][i];
						    		}
				 	    			
			 	    			}
				 	    			
				 	    			
				 	    			
				 	    			
			 	    			dc=Math.sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
						    	uhpy=1d/dc;
						    	sxn=pos[0]/dc;
			    				syn=pos[1]/dc;
			    				szn=pos[2]/dc;
			    				
			    				for(j=0;j<3;j++) {
				    				dotp=sxn*x[j][0]+syn*x[j][1]+szn*x[j][2];	
				    				
				    				vectmpx[0]=sxn*dotp;
				    				vectmpx[1]=syn*dotp;
				    				vectmpx[2]=szn*dotp;
				    				
				    				x[j][0]-=2d*vectmpx[0];
				    				x[j][1]-=2d*vectmpx[1];
				    				x[j][2]-=2d*vectmpx[2];
			    				}
			    				
			    				sxn=-pos[0]/dc;
			    				syn=-pos[1]/dc;
			    				szn=-pos[2]/dc;
			    				
			    				dotp=sxn*x[1][0]+syn*x[1][1]+szn*x[1][2];
			    				leangle=Math.acos(dotp);
			    				
			    				theta=Math.PI/2d-leangle;
			    					
				    				dist1=arsinh(Math.tan(theta));
				    				ang=Math.atan(Math.sinh(distrem-dist1));
				    				
				    				lerayon=uhpy/Math.sin(Math.PI/2d+theta);
				    				
				    				lex1=lerayon*Math.cos(Math.PI/2d+theta);
				    				
				    				lex2=lerayon*Math.cos(Math.PI/2d-ang);
				    				ley2=lerayon*Math.sin(Math.PI/2d-ang);
				    				
				    				sxn=-sxn;
				    				syn=-syn;
				    				szn=-szn;
				    				dotp=-dotp;
				    				leangle=Math.acos(dotp);
				    				
				    				vecperp[0]=x[1][0]-dotp*sxn;
				 	    			vecperp[1]=x[1][1]-dotp*syn;
				 	    			vecperp[2]=x[1][2]-dotp*szn;
				 	    				
				 	    			vecperpn = Math.sqrt(vecperp[0]*vecperp[0]+vecperp[1]*vecperp[1]+vecperp[2]*vecperp[2]);
				 	    				
				 	    			vecperp[0]/=vecperpn;
				 	    			vecperp[1]/=vecperpn;
				 	    			vecperp[2]/=vecperpn;

				 	    			
					 	    				exitlgt=lex2-lex1;
						 	    			exitangle=exitlgt%bhsize;
						 	    			exitangle=(2d*Math.PI/bhsize)*exitangle;
						 	    			newang=Math.PI/2d+ang;
						 	    			
						 	    			pos[0]=Math.cos(exitangle)*sxn+Math.sin(exitangle)*vecperp[0];
						 	    			pos[1]=Math.cos(exitangle)*syn+Math.sin(exitangle)*vecperp[1];
						 	    			pos[2]=Math.cos(exitangle)*szn+Math.sin(exitangle)*vecperp[2];
						 	    			
						 	    			pos[0]*=1/ley2;
						 	    			pos[1]*=1/ley2;
						 	    			pos[2]*=1/ley2;	
					 	    				
					 	    			newv[0]=Math.cos(Math.PI+exitangle-newang)*sxn+Math.sin(Math.PI+exitangle-newang)*vecperp[0];
					 	    			newv[1]=Math.cos(Math.PI+exitangle-newang)*syn+Math.sin(Math.PI+exitangle-newang)*vecperp[1];
					 	    			newv[2]=Math.cos(Math.PI+exitangle-newang)*szn+Math.sin(Math.PI+exitangle-newang)*vecperp[2];
					 	    			
					 	    			
					 	    			for(j=0;j<2;j++) 
					 	    			{
						 	    			dotp1=newv[0]*x[1][0]+newv[1]*x[1][1]+newv[2]*x[1][2];
						 	    			dotp2=newv[0]*x[others[1][j]][0]+newv[1]*x[others[1][j]][1]+newv[2]*x[others[1][j]][2];
						 	    			
						 	    			proj[0]=x[1][0]*dotp1+x[others[1][j]][0]*dotp2;
						 	    			proj[1]=x[1][1]*dotp1+x[others[1][j]][1]*dotp2;
						 	    			proj[2]=x[1][2]*dotp1+x[others[1][j]][2]*dotp2;
						 	    			
						 	    			projn = Math.sqrt(proj[0]*proj[0]+proj[1]*proj[1]+proj[2]*proj[2]);
					 	    				
						 	    			proj[0]/=projn;
						 	    			proj[1]/=projn;
						 	    			proj[2]/=projn;
						 	    			
						 	    			tmpl=proj[0]*x[1][0]+proj[1]*x[1][1]+proj[2]*x[1][2];
						 	    			if(Math.abs(tmpl)>1d) ang1=0;
						 	    			else
						 	    			{
							 	    			ang1=Math.acos(tmpl);
							 	    			ang1*=Math.signum(proj[0]*x[others[1][j]][0]+proj[1]*x[others[1][j]][1]+proj[2]*x[others[1][j]][2]);
						 	    			}
						 	    			for(i=0;i<3;i++)
								    		{
								    			newx[1][i]=x[1][i]*Math.cos(ang1)+Math.sin(ang1)*x[others[1][j]][i];
								    			newx[others[1][j]][i]=x[others[1][j]][i]*Math.cos(ang1)-Math.sin(ang1)*x[1][i];
								    			x[1][i]=newx[1][i];
								    			x[others[1][j]][i]=newx[others[1][j]][i];
								    		}
						 	    			
					 	    			}
			 	    			}
			 	    			x[1][0]*=-1d;
			    				x[1][1]*=-1d;
			    				x[1][2]*=-1d;
			 	    			
		    				}
		    		}
		    		if(holdd) 
		    		{
		    			dc=Math.sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
				    	uhpy=1d/dc;
				    	
		    			if(dc>1d) {
		    			 for(i=0;i<n;i++) pos[i]+=x[1][i]/10.0d;
		    			}
		    			else
		    			{
		    				
		    				sxn=-pos[0]/dc;
		    				syn=-pos[1]/dc;
		    				szn=-pos[2]/dc;
		    				
		    				dotp=sxn*x[1][0]+syn*x[1][1]+szn*x[1][2];
		    				leangle=Math.acos(dotp);
		    				
		    				theta=Math.PI/2d-leangle;
		    					
			    				dist1=arsinh(Math.tan(theta));
			    				ang=Math.atan(Math.sinh(0.1d-dist1));
			    				
			    				lerayon=uhpy/Math.sin(Math.PI/2d+theta);
			    				
			    				lex1=lerayon*Math.cos(Math.PI/2d+theta);
			    				
			    				lex2=lerayon*Math.cos(Math.PI/2d-ang);
			    				ley2=lerayon*Math.sin(Math.PI/2d-ang);
			    				
			    				sxn=-sxn;
			    				syn=-syn;
			    				szn=-szn;
			    				dotp=-dotp;
			    				leangle=Math.acos(dotp);
			    				
			    				vecperp[0]=x[1][0]-dotp*sxn;
			 	    			vecperp[1]=x[1][1]-dotp*syn;
			 	    			vecperp[2]=x[1][2]-dotp*szn;
			 	    				
			 	    			vecperpn = Math.sqrt(vecperp[0]*vecperp[0]+vecperp[1]*vecperp[1]+vecperp[2]*vecperp[2]);
			 	    				
			 	    			vecperp[0]/=vecperpn;
			 	    			vecperp[1]/=vecperpn;
			 	    			vecperp[2]/=vecperpn;
			 	    			
			 	    			
			 	    			
			 	    			if(ley2<whsize)
			 	    			{
			 	    			
				 	    				exitlgt=lex2-lex1;
					 	    			exitangle=exitlgt%bhsize;
					 	    			exitangle=(2d*Math.PI/bhsize)*exitangle;
					 	    			newang=Math.PI/2d+ang;
					 	    			
					 	    			pos[0]=Math.cos(exitangle)*sxn+Math.sin(exitangle)*vecperp[0];
					 	    			pos[1]=Math.cos(exitangle)*syn+Math.sin(exitangle)*vecperp[1];
					 	    			pos[2]=Math.cos(exitangle)*szn+Math.sin(exitangle)*vecperp[2];
					 	    			
					 	    			pos[0]*=1/ley2;
					 	    			pos[1]*=1/ley2;
					 	    			pos[2]*=1/ley2;	
				 	    				
				 	    			newv[0]=Math.cos(Math.PI+exitangle-newang)*sxn+Math.sin(Math.PI+exitangle-newang)*vecperp[0];
				 	    			newv[1]=Math.cos(Math.PI+exitangle-newang)*syn+Math.sin(Math.PI+exitangle-newang)*vecperp[1];
				 	    			newv[2]=Math.cos(Math.PI+exitangle-newang)*szn+Math.sin(Math.PI+exitangle-newang)*vecperp[2];
				 	    			
				 	    			
				 	    			for(j=0;j<2;j++) 
				 	    			{
					 	    			dotp1=newv[0]*x[1][0]+newv[1]*x[1][1]+newv[2]*x[1][2];
					 	    			dotp2=newv[0]*x[others[1][j]][0]+newv[1]*x[others[1][j]][1]+newv[2]*x[others[1][j]][2];
					 	    			
					 	    			proj[0]=x[1][0]*dotp1+x[others[1][j]][0]*dotp2;
					 	    			proj[1]=x[1][1]*dotp1+x[others[1][j]][1]*dotp2;
					 	    			proj[2]=x[1][2]*dotp1+x[others[1][j]][2]*dotp2;
					 	    			
					 	    			projn = Math.sqrt(proj[0]*proj[0]+proj[1]*proj[1]+proj[2]*proj[2]);
				 	    				
					 	    			proj[0]/=projn;
					 	    			proj[1]/=projn;
					 	    			proj[2]/=projn;
					 	    			
					 	    			tmpl=proj[0]*x[1][0]+proj[1]*x[1][1]+proj[2]*x[1][2];
					 	    			if(Math.abs(tmpl)>1d) ang1=0;
					 	    			else
					 	    			{
						 	    			ang1=Math.acos(tmpl);
						 	    			ang1*=Math.signum(proj[0]*x[others[1][j]][0]+proj[1]*x[others[1][j]][1]+proj[2]*x[others[1][j]][2]);
					 	    			}
					 	    			for(i=0;i<3;i++)
							    		{
							    			newx[1][i]=x[1][i]*Math.cos(ang1)+Math.sin(ang1)*x[others[1][j]][i];
							    			newx[others[1][j]][i]=x[others[1][j]][i]*Math.cos(ang1)-Math.sin(ang1)*x[1][i];
							    			x[1][i]=newx[1][i];
							    			x[others[1][j]][i]=newx[others[1][j]][i];
							    		}
					 	    			
				 	    			}
			 	    			}
			 	    			else
			 	    			{
			 	    				
			 	    				flipg=!flipg;
	
			 	    				lex2=-Math.sqrt(lerayon*lerayon-whsize*whsize);
			 	    				
			 	    				angx1=Math.PI/2d - Math.asin(whsize/lerayon);
				 	    			dist2=arsinh(Math.tan(angx1));
				 	    			distrem=0.1d-dist1+dist2;
				 	    			ang=Math.atan(Math.sinh(0.1d-distrem-dist1));
				 	    			
				 	    			exitlgt=lex2-lex1;
				 	    			exitangle=exitlgt%bhsize;
				 	    			exitangle=(2d*Math.PI/bhsize)*exitangle;
				 	    			newang=Math.PI/2d+ang;
				 	    			
				 	    			pos[0]=Math.cos(exitangle)*sxn+Math.sin(exitangle)*vecperp[0];
				 	    			pos[1]=Math.cos(exitangle)*syn+Math.sin(exitangle)*vecperp[1];
				 	    			pos[2]=Math.cos(exitangle)*szn+Math.sin(exitangle)*vecperp[2];
				 	    			
				 	    			pos[0]*=1/whsize;
				 	    			pos[1]*=1/whsize;
				 	    			pos[2]*=1/whsize;	
			 	    				
			 	    			newv[0]=Math.cos(Math.PI+exitangle-newang)*sxn+Math.sin(Math.PI+exitangle-newang)*vecperp[0];
			 	    			newv[1]=Math.cos(Math.PI+exitangle-newang)*syn+Math.sin(Math.PI+exitangle-newang)*vecperp[1];
			 	    			newv[2]=Math.cos(Math.PI+exitangle-newang)*szn+Math.sin(Math.PI+exitangle-newang)*vecperp[2];
			 	    			
			 	    			
			 	    			for(j=0;j<2;j++) 
			 	    			{
				 	    			dotp1=newv[0]*x[1][0]+newv[1]*x[1][1]+newv[2]*x[1][2];
				 	    			dotp2=newv[0]*x[others[1][j]][0]+newv[1]*x[others[1][j]][1]+newv[2]*x[others[1][j]][2];
				 	    			
				 	    			proj[0]=x[1][0]*dotp1+x[others[1][j]][0]*dotp2;
				 	    			proj[1]=x[1][1]*dotp1+x[others[1][j]][1]*dotp2;
				 	    			proj[2]=x[1][2]*dotp1+x[others[1][j]][2]*dotp2;
				 	    			
				 	    			projn = Math.sqrt(proj[0]*proj[0]+proj[1]*proj[1]+proj[2]*proj[2]);
			 	    				
				 	    			proj[0]/=projn;
				 	    			proj[1]/=projn;
				 	    			proj[2]/=projn;
				 	    			
				 	    			tmpl=proj[0]*x[1][0]+proj[1]*x[1][1]+proj[2]*x[1][2];
				 	    			if(Math.abs(tmpl)>1d) ang1=0;
				 	    			else
				 	    			{
					 	    			ang1=Math.acos(tmpl);
					 	    			ang1*=Math.signum(proj[0]*x[others[1][j]][0]+proj[1]*x[others[1][j]][1]+proj[2]*x[others[1][j]][2]);
				 	    			}
				 	    			for(i=0;i<3;i++)
						    		{
						    			newx[1][i]=x[1][i]*Math.cos(ang1)+Math.sin(ang1)*x[others[1][j]][i];
						    			newx[others[1][j]][i]=x[others[1][j]][i]*Math.cos(ang1)-Math.sin(ang1)*x[1][i];
						    			x[1][i]=newx[1][i];
						    			x[others[1][j]][i]=newx[others[1][j]][i];
						    		}
				 	    			
			 	    			}
				 	    			
				 	    			
				 	    			
				 	    			
			 	    			dc=Math.sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
						    	uhpy=1d/dc;
						    	sxn=pos[0]/dc;
			    				syn=pos[1]/dc;
			    				szn=pos[2]/dc;
			    				
			    				for(j=0;j<3;j++) {
				    				dotp=sxn*x[j][0]+syn*x[j][1]+szn*x[j][2];	
				    				
				    				vectmpx[0]=sxn*dotp;
				    				vectmpx[1]=syn*dotp;
				    				vectmpx[2]=szn*dotp;
				    				
				    				x[j][0]-=2d*vectmpx[0];
				    				x[j][1]-=2d*vectmpx[1];
				    				x[j][2]-=2d*vectmpx[2];
			    				}
			    				
			    				sxn=-pos[0]/dc;
			    				syn=-pos[1]/dc;
			    				szn=-pos[2]/dc;
			    				
			    				dotp=sxn*x[1][0]+syn*x[1][1]+szn*x[1][2];
			    				leangle=Math.acos(dotp);
			    				
			    				theta=Math.PI/2d-leangle;
			    					
				    				dist1=arsinh(Math.tan(theta));
				    				ang=Math.atan(Math.sinh(distrem-dist1));
				    				
				    				lerayon=uhpy/Math.sin(Math.PI/2d+theta);
				    				
				    				lex1=lerayon*Math.cos(Math.PI/2d+theta);
				    				
				    				lex2=lerayon*Math.cos(Math.PI/2d-ang);
				    				ley2=lerayon*Math.sin(Math.PI/2d-ang);
				    				
				    				sxn=-sxn;
				    				syn=-syn;
				    				szn=-szn;
				    				dotp=-dotp;
				    				leangle=Math.acos(dotp);
				    				
				    				vecperp[0]=x[1][0]-dotp*sxn;
				 	    			vecperp[1]=x[1][1]-dotp*syn;
				 	    			vecperp[2]=x[1][2]-dotp*szn;
				 	    				
				 	    			vecperpn = Math.sqrt(vecperp[0]*vecperp[0]+vecperp[1]*vecperp[1]+vecperp[2]*vecperp[2]);
				 	    				
				 	    			vecperp[0]/=vecperpn;
				 	    			vecperp[1]/=vecperpn;
				 	    			vecperp[2]/=vecperpn;

				 	    			
					 	    				exitlgt=lex2-lex1;
						 	    			exitangle=exitlgt%bhsize;
						 	    			exitangle=(2d*Math.PI/bhsize)*exitangle;
						 	    			newang=Math.PI/2d+ang;
						 	    			
						 	    			pos[0]=Math.cos(exitangle)*sxn+Math.sin(exitangle)*vecperp[0];
						 	    			pos[1]=Math.cos(exitangle)*syn+Math.sin(exitangle)*vecperp[1];
						 	    			pos[2]=Math.cos(exitangle)*szn+Math.sin(exitangle)*vecperp[2];
						 	    			
						 	    			pos[0]*=1/ley2;
						 	    			pos[1]*=1/ley2;
						 	    			pos[2]*=1/ley2;	
					 	    				
					 	    			newv[0]=Math.cos(Math.PI+exitangle-newang)*sxn+Math.sin(Math.PI+exitangle-newang)*vecperp[0];
					 	    			newv[1]=Math.cos(Math.PI+exitangle-newang)*syn+Math.sin(Math.PI+exitangle-newang)*vecperp[1];
					 	    			newv[2]=Math.cos(Math.PI+exitangle-newang)*szn+Math.sin(Math.PI+exitangle-newang)*vecperp[2];
					 	    			
					 	    			
					 	    			for(j=0;j<2;j++) 
					 	    			{
						 	    			dotp1=newv[0]*x[1][0]+newv[1]*x[1][1]+newv[2]*x[1][2];
						 	    			dotp2=newv[0]*x[others[1][j]][0]+newv[1]*x[others[1][j]][1]+newv[2]*x[others[1][j]][2];
						 	    			
						 	    			proj[0]=x[1][0]*dotp1+x[others[1][j]][0]*dotp2;
						 	    			proj[1]=x[1][1]*dotp1+x[others[1][j]][1]*dotp2;
						 	    			proj[2]=x[1][2]*dotp1+x[others[1][j]][2]*dotp2;
						 	    			
						 	    			projn = Math.sqrt(proj[0]*proj[0]+proj[1]*proj[1]+proj[2]*proj[2]);
					 	    				
						 	    			proj[0]/=projn;
						 	    			proj[1]/=projn;
						 	    			proj[2]/=projn;
						 	    			
						 	    			tmpl=proj[0]*x[1][0]+proj[1]*x[1][1]+proj[2]*x[1][2];
						 	    			if(Math.abs(tmpl)>1d) ang1=0;
						 	    			else
						 	    			{
							 	    			ang1=Math.acos(tmpl);
							 	    			ang1*=Math.signum(proj[0]*x[others[1][j]][0]+proj[1]*x[others[1][j]][1]+proj[2]*x[others[1][j]][2]);
						 	    			}
						 	    			for(i=0;i<3;i++)
								    		{
								    			newx[1][i]=x[1][i]*Math.cos(ang1)+Math.sin(ang1)*x[others[1][j]][i];
								    			newx[others[1][j]][i]=x[others[1][j]][i]*Math.cos(ang1)-Math.sin(ang1)*x[1][i];
								    			x[1][i]=newx[1][i];
								    			x[others[1][j]][i]=newx[others[1][j]][i];
								    		}
						 	    			
					 	    			}
			 	    			}
			 	    			
			 	    			
		    				}
		    		}
		    		
		    		
		    		
		    		mousx=MouseInfo.getPointerInfo().getLocation().x-centralx;
		    		mousy=MouseInfo.getPointerInfo().getLocation().y-centraly;
		    		
		    		if(mousx>0)
		    		{
		    			
		    			
			    			for(i=0;i<n;i++)
			    			{
			    				newx[0][i]=x[0][i]*rcos+rsin*x[1][i];
			    				newx[1][i]=x[1][i]*rcos-rsin*x[0][i];
			    				x[0][i]=newx[0][i];
			    				x[1][i]=newx[1][i];
			    			}
		    			
		    		}
		    		else if(mousx<0)
		    		{
		    			
		    			
			    			for(i=0;i<n;i++)
			    			{
			    				newx[0][i]=x[0][i]*rcos-rsin*x[1][i];
			    				newx[1][i]=x[1][i]*rcos+rsin*x[0][i];
			    				x[0][i]=newx[0][i];
			    				x[1][i]=newx[1][i];
			    			}
		    			
		    		}
		    		
		    		if(mousy<0)
		    		{
		    			
		    		
			    			for(i=0;i<n;i++)
			    			{
			    				newx[0][i]=x[0][i]*rcos+rsin*x[2][i];
			    				newx[2][i]=x[2][i]*rcos-rsin*x[0][i];
			    				x[0][i]=newx[0][i];
			    				x[2][i]=newx[2][i];
			    			}
		    			
		    		}
		    		else if(mousy>0)
		    		{
		    		
			    			for(i=0;i<n;i++)
			    			{
			    				newx[0][i]=x[0][i]*rcos-rsin*x[2][i];
			    				newx[2][i]=x[2][i]*rcos+rsin*x[0][i];
			    				x[0][i]=newx[0][i];
			    				x[2][i]=newx[2][i];
			    			}
		    			
		    		}
		    		
			    	currentpix=0;
			    	for(i=0;i<n;i++)
			    	{
			    		vec[i]=dist*x[0][i]+multy*x[1][i]+multz*x[2][i];
			    		addy[i]=sqsz*x[1][i];
			    		addz[i]=msqsz*x[2][i];
			    	}
			    	
			    	
			    	for(i=0;i<n;i++)
			    	{
			    		vec[i]=dist*x[0][i]+multy*x[1][i]+multz*x[2][i];
			    		addy[i]=sqsz*x[1][i];
			    		addz[i]=msqsz*x[2][i];
			    	}
			    	
			    	dc=Math.sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
			    	uhpy=1d/dc;
				    	
			 	    for(ph=0;ph<height;ph++)
			 	    {
			 	    	for(i=0;i<n;i++) vectmp[i]=vec[i];
			 	    	for(pw=0;pw<width;pw++)
			 	    	{
		
			 	    		for(i=0;i<n;i++)
			 				{
			 					vecn[i]=vec[i]*vecl[ph][pw];
			 				}
			 	    		
			 	    		flip=false;
			 				 	    		
			 				 	    		if(dc<1d)
			 				 	    		{
			 				 	    			sx=pos[0]/dc;
			 				 	    			sy=pos[1]/dc;
			 				 	    			sz=pos[2]/dc;
			 				 	    			
			 				 	    			dotp=sx*vecn[0]+sy*vecn[1]+sz*vecn[2];
			 			 	    				leangle=Math.acos(dotp);
			 			 	    					
			 			 	    				vecperp[0]=vecn[0]-dotp*sx;
			 			 	    				vecperp[1]=vecn[1]-dotp*sy;
			 			 	    				vecperp[2]=vecn[2]-dotp*sz;
			 			 	    				
			 			 	    				vecperpn = Math.sqrt(vecperp[0]*vecperp[0]+vecperp[1]*vecperp[1]+vecperp[2]*vecperp[2]);
			 			 	    				
			 			 	    				vecperp[0]/=vecperpn;
			 			 	    				vecperp[1]/=vecperpn;
			 			 	    				vecperp[2]/=vecperpn;
			 			 	    				
			 			 	    				lecos=Math.cos(leangle);
			 			 	    				lesign=Math.signum(lecos);
			 			 	    				
			 			 	    				lecossq=lecos*lecos;
			 			 	    				
			 			 	    				
			 			 	    				lerayon=uhpy*Math.sqrt(1d/(1-lecossq));
			 			 	    				
			 			 	    				
				 			 	    				lex1=Math.sqrt(lerayon*lerayon - uhpy*uhpy)*lesign;
				 			 	    				lex2=Math.sqrt(lerayon*lerayon-1);
				 			 	    				
				 			 	    			
				 			 	    				
				 			 	    				if(lesign<0 && lerayon>whsize)
				 			 	    				{
				 			 	    					flip=true;
				 			 	    					
				 			 	    					lexw=Math.sqrt(lerayon*lerayon-whsize*whsize);
				 			 	    					
				 			 	    					exitlgt=-lex1+lex2-2d*lexw;
					 			 	    				exitangle=exitlgt%bhsize;
					 			 	    				exitangle=(2d*Math.PI/bhsize)*exitangle;
					 			 	    				
					 			 	    				exitangle2=Math.acos((1d/uhpy)*Math.sqrt(uhpy*uhpy-1d+lecossq));
					 			 	    				
					 			 	    				
				 			 	    					pos2[0]=Math.cos(exitangle)*sx+Math.sin(exitangle)*vecperp[0];
					 			 	    				pos2[1]=Math.cos(exitangle)*sy+Math.sin(exitangle)*vecperp[1];
					 			 	    				pos2[2]=Math.cos(exitangle)*sz+Math.sin(exitangle)*vecperp[2];
					 			 	    				
					 			 	    				vecn2[0]=Math.cos(exitangle+exitangle2)*sx+Math.sin(exitangle2+exitangle)*vecperp[0];
					 			 	    				vecn2[1]=Math.cos(exitangle+exitangle2)*sy+Math.sin(exitangle2+exitangle)*vecperp[1];
					 			 	    				vecn2[2]=Math.cos(exitangle+exitangle2)*sz+Math.sin(exitangle2+exitangle)*vecperp[2];
				 			 	    				}
				 			 	    				else
				 			 	    				{
				 			 	    					exitlgt=lex2-lex1;
					 			 	    				exitangle=exitlgt%bhsize;
					 			 	    				exitangle=(2d*Math.PI/bhsize)*exitangle;
					 			 	    				
					 			 	    				exitangle2=Math.acos((1d/uhpy)*Math.sqrt(uhpy*uhpy-1d+lecossq));
					 			 	    				
					 			 	    				pos2[0]=Math.cos(exitangle)*sx+Math.sin(exitangle)*vecperp[0];
					 			 	    				pos2[1]=Math.cos(exitangle)*sy+Math.sin(exitangle)*vecperp[1];
					 			 	    				pos2[2]=Math.cos(exitangle)*sz+Math.sin(exitangle)*vecperp[2];
		
					 			 	    				
					 			 	    				vecn2[0]=Math.cos(exitangle+exitangle2)*sx+Math.sin(exitangle2+exitangle)*vecperp[0];
					 			 	    				vecn2[1]=Math.cos(exitangle+exitangle2)*sy+Math.sin(exitangle2+exitangle)*vecperp[1];
					 			 	    				vecn2[2]=Math.cos(exitangle+exitangle2)*sz+Math.sin(exitangle2+exitangle)*vecperp[2];
				 			 	    				}
				 			 	    				//tf=2d*Math.log((Math.sqrt((lex2-lex1)*(lex2-lex1)+(1-uhpy)*(1-uhpy))+Math.sqrt((lex2-lex1)*(lex2-lex1)+(1+uhpy)*(1+uhpy)))/(2d*Math.sqrt(uhpy)));
				 			 	    				tcont=0;
				 			 	    				
					 			 	    				
			 				 	    		}
			 				 	    		else
			 				 	    		{

			 				 	    			qa=vecn[0]*vecn[0]+vecn[1]*vecn[1]+vecn[2]*vecn[2];
			 				 	    			qb=2d*(vecn[0]*pos[0]+vecn[1]*pos[1]+vecn[2]*pos[2]);
			 				 	    			qc=pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]-1d;
			 				 	    			
			 				 	    			discr=qb*qb-4d*qa*qc;
			 				 	    			
			 				 	    			if(discr<=0) {
			 					 	    			pos2[0]=pos[0];
			 					 	    			pos2[1]=pos[1];
			 					 	    			pos2[2]=pos[2];
			 					 	    			
			 					 	    			vecn2[0]=vecn[0];
			 					 	    			vecn2[1]=vecn[1];
			 					 	    			vecn2[2]=vecn[2];
			 					 	    			
			 					 	    			tcont=0;
			 					 	    			//tf=0;
			 				 	    			}
			 				 	    			else
			 				 	    			{
			 				 	    				
			 				 	    				tcont=((-1d)*qb-Math.sqrt(discr))/(2d*qa);
			 				 	    				
			 				 	    				if(tcont>0) {
			 				 	    				sx=vecn[0]*tcont+pos[0];
			 				 	    				sy=vecn[1]*tcont+pos[1];
			 				 	    				sz=vecn[2]*tcont+pos[2];
			 				 	    				
			 				 	    				dotp=sx*vecn[0]+sy*vecn[1]+sz*vecn[2];
			 				 	    				leangle=Math.acos(dotp);
			 				 	    					
			 				 	    				vecperp[0]=vecn[0]-dotp*sx;
			 				 	    				vecperp[1]=vecn[1]-dotp*sy;
			 				 	    				vecperp[2]=vecn[2]-dotp*sz;
			 				 	    				
			 				 	    				vecperpn = Math.sqrt(vecperp[0]*vecperp[0]+vecperp[1]*vecperp[1]+vecperp[2]*vecperp[2]);
			 				 	    				
			 				 	    				vecperp[0]/=vecperpn;
			 				 	    				vecperp[1]/=vecperpn;
			 				 	    				vecperp[2]/=vecperpn;
			 				 	    				
			 				 	    				lecos=Math.cos(leangle);
			 				 	    				lecossq=lecos*lecos;
			 				 	    			
			 				 	    				lerayon=Math.sqrt(1d/(1-lecossq));
			 				 	    				
			 				 	    				
			 				 	    				
			 				 	    				
			 				 	    				if(lerayon>whsize)
			 				 	    				{
			 				 	    					
			 				 	    					flip=true;
			 				 	    					
			 				 	    					
			 				 	    					lexw=Math.sqrt((1d/(1-lecossq))-whsize*whsize);
			 				 	    					
			 				 	    					exitlgt=-2d*(lexw+lerayon*lecos);
				 				 	    				exitangle=exitlgt%bhsize;
				 			
				 				 	    				exitangle=(2d*Math.PI/bhsize)*exitangle;
			 				 	    					
				 				 	    				pos2[0]=Math.cos(exitangle)*sx+Math.sin(exitangle)*vecperp[0];
				 				 	    				pos2[1]=Math.cos(exitangle)*sy+Math.sin(exitangle)*vecperp[1];
				 				 	    				pos2[2]=Math.cos(exitangle)*sz+Math.sin(exitangle)*vecperp[2];
	
				 				 	    				
				 				 	    				vecn2[0]=Math.cos(Math.PI+exitangle-leangle)*sx+Math.sin(Math.PI-leangle+exitangle)*vecperp[0];
				 				 	    				vecn2[1]=Math.cos(Math.PI+exitangle-leangle)*sy+Math.sin(Math.PI-leangle+exitangle)*vecperp[1];
				 				 	    				vecn2[2]=Math.cos(Math.PI+exitangle-leangle)*sz+Math.sin(Math.PI-leangle+exitangle)*vecperp[2];
			 				 	    				}
			 				 	    				else
			 				 	    				{
				 				 	    				exitlgt=-2d*lerayon*lecos;
				 				 	    				exitangle=exitlgt%bhsize;
				 			
				 				 	    				exitangle=(2d*Math.PI/bhsize)*exitangle;
				 				 	    				
	
				 				 	    				pos2[0]=Math.cos(exitangle)*sx+Math.sin(exitangle)*vecperp[0];
				 				 	    				pos2[1]=Math.cos(exitangle)*sy+Math.sin(exitangle)*vecperp[1];
				 				 	    				pos2[2]=Math.cos(exitangle)*sz+Math.sin(exitangle)*vecperp[2];
	
				 				 	    				
				 				 	    				vecn2[0]=Math.cos(Math.PI+exitangle-leangle)*sx+Math.sin(Math.PI-leangle+exitangle)*vecperp[0];
				 				 	    				vecn2[1]=Math.cos(Math.PI+exitangle-leangle)*sy+Math.sin(Math.PI-leangle+exitangle)*vecperp[1];
				 				 	    				vecn2[2]=Math.cos(Math.PI+exitangle-leangle)*sz+Math.sin(Math.PI-leangle+exitangle)*vecperp[2];
			 				 	    				}
			 				 	    			    //tf=2d*Math.log((exitlgt+Math.sqrt(exitlgt*exitlgt+4))/2d);
			 				 	    				}
			 				 	    				else {
			 				 	    					pos2[0]=pos[0];
				 					 	    			pos2[1]=pos[1];
				 					 	    			pos2[2]=pos[2];
				 					 	    			
				 					 	    			vecn2[0]=vecn[0];
				 					 	    			vecn2[1]=vecn[1];
				 					 	    			vecn2[2]=vecn[2];
				 					 	    			
				 					 	    			tcont=0;
				 					 	    			//tf=0;
			 				 	    				}

			 				 	    			}
			 				 	    		}
			 				 	    		
			 				 	    		
			 				 	    		tmin=(Math.signum(vecn2[0])*roomsize-pos2[0])/vecn2[0];
			 			 	    			tmincoord=0;
			 			 	    			
			 			 	    			tsol=(Math.signum(vecn2[1])*roomsize-pos2[1])/vecn2[1];
			 			 	    			if(tsol<tmin)
			 			 	    			{
			 			 	    				tmin=tsol;
			 			 	    				tmincoord=1;
			 			 	    			}
			 			 	    			
			 			 	    			tsol=(Math.signum(vecn2[2])*roomsize-pos2[2])/vecn2[2];
			 			 	    			if(tsol<tmin)
			 			 	    			{
			 			 	    				tmin=tsol;
			 			 	    				tmincoord=2;
			 			 	    			}
			 			 	    			
			 			 	    			setcoord=0;
			 			 	    			for(i=0;i<3;i++)
			 			 	    			{
			 			 	    				if(i!=tmincoord)
			 			 	    				{
			 			 	    					coll[setcoord]=pos2[i]+tmin*vecn2[i];
			 			 	    					setcoord++;
			 			 	    				}
			 			 	    			}
			 			 	    			
			 			 	    			checker=((int)Math.floor(coll[0]*schecker))%2;
			 			 	    			checker+=((int)Math.floor(coll[1]*schecker))%2;
			 			 	    			if(checker<0) checker+=2;
			 			 	    			checker%=2;
			 			 	    			
			 			 	    			
			 			 	    			
			 			 	    			
			 			 	    			
			 			 	    			  //ttot=tf+tcont+tmin;
			 			 	    			  //alpha=((-1d)/seedist)*ttot+1;
			 			 	    			  //if(alpha<0) alpha=0d;
			 			 	    	
			 			 	    			alpha=255d;
			 			 	    			
			 			 	    			if(flip^flipg)
			 			 	    			{
			 			 	    				if(tmincoord==0)
		 			 	    					{
		 			 	    						ctmp[0]=0;
		 			 	    						ctmp[1]=0;
		 			 	    						ctmp[2]=255;
		 			 	    					}
			 			 	    				else if(tmincoord==1)
			 			 	    				{
			 			 	    					ctmp[0]=255;
			 			 	    					ctmp[1]=0;
			 			 	    					ctmp[2]=0;
			 			 	    				}
			 			 	    				else if(tmincoord==2)
			 			 	    				{
			 			 	    					ctmp[0]=0;
			 			 	    					ctmp[1]=255;
			 			 	    					ctmp[2]=255;
			 			 	    				}
			 			 	    				
			 			 	    				
			 			 	    				
			 			 	    				if(checker==0) {
			 			 	    					
			 			 	    					
			 			 	    					
				 			 	    				pixels[currentpix]=(byte)alpha;
				 					 	    		currentpix++;
				 					 	    		pixels[currentpix]=(byte)0;
				 					 	    		currentpix++;
				 					 	    		pixels[currentpix]=(byte)0;
				 					 	    		currentpix++;
				 					 	    		pixels[currentpix]=(byte)0;
				 					 	    		currentpix++;
				 			 	    			}
				 			 	    			else
				 			 	    			{
				 			 	    				pixels[currentpix]=(byte)alpha;
				 					 	    		currentpix++;
				 			 	    				pixels[currentpix]=(byte)ctmp[0];
				 					 	    		currentpix++;
				 					 	    		pixels[currentpix]=(byte)ctmp[1];
				 					 	    		currentpix++;
				 					 	    		pixels[currentpix]=(byte)ctmp[2];
				 					 	    		currentpix++;
				 			 	    			}
			 			 	    			}
			 			 	    			else
			 			 	    			{
			 			 	    				if(tmincoord==2) {
					 			 	    			ctmp[0]=255;
					 			 	    			ctmp[1]=255;
					 			 	    			ctmp[2]=255; }
					 			 	    			else if(tmincoord==0)
					 			 	    			{
					 			 	    				if(vecn2[0]<0) ctmp=rnbw((1d/8d)-coll[0]/(8d*roomsize));
					 			 	    				else ctmp=rnbw((1d/2d)+(1d/8d)+coll[0]/(8d*roomsize)); 
					 			 	    			}
					 			 	    			else
					 			 	    			{
					 			 	    				if(vecn2[1]<0) ctmp=rnbw((1d/4d)+(1d/8d)+coll[0]/(8d*roomsize));
					 			 	    				else ctmp=rnbw((3d/4d)+(1d/8d)-coll[0]/(8d*roomsize));
					 			 	    			}
			 			 	    				
				 			 	    			if(checker==0) {
				 			 	    				pixels[currentpix]=(byte)alpha;
				 					 	    		currentpix++;
				 					 	    		pixels[currentpix]=(byte)0;
				 					 	    		currentpix++;
				 					 	    		pixels[currentpix]=(byte)0;
				 					 	    		currentpix++;
				 					 	    		pixels[currentpix]=(byte)0;
				 					 	    		currentpix++;
				 			 	    			}
				 			 	    			else
				 			 	    			{
				 			 	    				pixels[currentpix]=(byte)alpha;
				 					 	    		currentpix++;
				 			 	    				pixels[currentpix]=(byte)ctmp[0];
				 					 	    		currentpix++;
				 					 	    		pixels[currentpix]=(byte)ctmp[1];
				 					 	    		currentpix++;
				 					 	    		pixels[currentpix]=(byte)ctmp[2];
				 					 	    		currentpix++;
				 			 	    			}
			 			 	    			}
			 
			 	    		
			 	    		for(i=0;i<n;i++) vec[i]+=addy[i];
			 	    	}
			 	    	
			 	    	for(i=0;i<n;i++) vec[i]=vectmp[i]+addz[i];
			 	    }
			 	    

			       label.setIcon(new ImageIcon(image));
			       frame.getContentPane().add(label,BorderLayout.CENTER);
			       robot.mouseMove(centralx,centraly);
	    	}
	    }
	}
	
	public static double vecprod(double[] v1, double[] v2)
	{
		double ret=0;
		for(int i=0;i<n;i++) ret+=v1[i]*v2[i];
		return ret;
	}
	
	
	public static int[] rnbw(double x)
	{
		int[] ret = new int[3];
		double tmp;
		
		tmp=x%(1d/6d);
		
		if(x<1d/6d)
		{
			ret[0]=255;
			ret[1]=(int) (1530d*tmp);
		}
		else if(x<1d/3d)
		{
			ret[1]=255;
			ret[0]=(int) (255d-1530d*tmp);
		}
		else if(x<0.5d)
		{
			ret[1]=255;
			ret[2]=(int) (1530d*tmp);
		}
		else if(x<2d/3d)
		{
			ret[2]=255;
			ret[1]=(int) (255d-1530d*tmp);
		}
		else if(x<5d/6d)
		{
			ret[2]=255;
			ret[0]=(int) (1530d*tmp);
		}
		else
		{
			ret[0]=255;
			ret[2]=(int) (255d-1530d*tmp);
		}
		
		return ret;
	}
	

	public static double veclgt(double[] v)
	{
		return (double)Math.sqrt(vecprod(v,v));
	}
	


	public static double arsinh(double x)
	{
		return Math.log(x+Math.sqrt(x*x+1));
	}
	
	@Override
	public void focusLost(FocusEvent e) {
        focus=false;
        frame.getContentPane().setCursor(Cursor.getDefaultCursor());
    }
	
	public whs() {
        addKeyListener(this);
        addFocusListener(this);
    }
	
	public void addNotify() {
        super.addNotify();
        requestFocus();
    }
	
	@Override
	public void keyPressed(KeyEvent e) { }
	@Override
	public void keyReleased(KeyEvent e) {
		c = e.getKeyChar();
		if(c==119) {holdw=false;}
		else if(c==97) {holda=false;}
		else if(c==115) {holds=false;}
		else if(c==100) {holdd=false;}
		else if(c==122) {holdz=false;}
		else if(c==120) {holdx=false;}
		else if(c==99) {holdc=false;}
		else if(c==118) {holdv=false;}
	}
	@Override
	public void keyTyped(KeyEvent e) {
		c = e.getKeyChar();
		if(c==27) { focus=false;frame.getContentPane().setCursor(Cursor.getDefaultCursor());}
		else if(c==119) {holdw=true;}
		else if(c==97) {holda=true;}
		else if(c==115) {holds=true;}
		else if(c==100) {holdd=true;}

		else if(c==112) {resettot=true;}
		else if(c==122) {holdz=true;}
		else if(c==120) {holdx=true;}
		else if(c==99) {holdc=true;}
		else if(c==118) {holdv=true;}
	}

	@Override
	public void focusGained(FocusEvent arg0) {}
}