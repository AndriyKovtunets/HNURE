package ua.nure;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Locale;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;


public class Main {
	static Random rand = new Random();
 

	static DecimalFormat df = new DecimalFormat("#.##");
	static DecimalFormatSymbols otherSymbols = new DecimalFormatSymbols(Locale.getDefault());

	final static float x1_min = 0.5f, x1_max = 1.1f, x2_min = 1.0f, x2_max = 4.6f;
	final static int n = 500, amount = 2;
	static int iterations=0;



	public static void main(String[] args) {
		//		 otherSymbols.setGroupingSeparator('.'); 
		otherSymbols.setDecimalSeparator('.');
		df = new DecimalFormat("#.##", otherSymbols);

		float instance[][] = new float[n][amount];

		//generate instance
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < amount; j++) {
				if (j == 0)
					instance[i][j] = Float.parseFloat(df.format( rand.nextFloat() * (x1_max - x1_min) + x1_min )) ;
				if (j == 1)
					instance[i][j] = Float.parseFloat(df.format( rand.nextFloat() * (x2_max - x2_min) + x2_min )) ;
			}
		}

		String[][] s_instance = get_StringBinaryInstance(n, amount, x1_min, x1_max, x2_min, x2_max);
		
		
		//  parent's choice
		long timet_start, timer_end;
		timet_start=System.nanoTime();
		
		ArrayList<String> crossedChromosomes = new ArrayList<String>(); 
		TreeMap<Float,Integer>  treeMap = new TreeMap<Float,Integer>(Collections.reverseOrder());

		TreeMap<Float,Float[]>  best = new TreeMap<Float,Float[]>();

		
		int count=(int) (n*0.9);
		int parent1 = 0,parent2=0;
		float[] best_rezult=new float[3];
		boolean mask[][]=new boolean[n][amount];
		String chromosome1 = null,chromosome2= null;
		//цыкл

		int ii;
		for(ii=0; ii<100; ii++){

			
		
		for (int i=0; i<count;i++){
			parent1=rand.nextInt(499);
//			mask[parent1][0]=true;
			chromosome1 = s_instance[parent1][0];
			//виб≥р parent2
			chromosome2=get_SimilarParent(chromosome1, s_instance, count);
			//схрещуванн€
			crossedChromosomes=get_crossedChromosomes(crossedChromosomes, chromosome1, chromosome2);
		}

		//велика хромосома)
		StringBuffer long_Chromosome = new StringBuffer();
		for(String s: crossedChromosomes){
			long_Chromosome.append(s);
		}

		
		
		// мутац≥€

		
		for(int i=0; i<long_Chromosome.length();i++){
			if(rand.nextDouble()<=0.01 )  {long_Chromosome.setCharAt(i, long_Chromosome.charAt(i)==0 ? '1':'0'); }
		}



		//розбиваэм велику хромосому
		crossedChromosomes.clear();
		for(int i=0; i<=long_Chromosome.length()-30; i+=30){

			crossedChromosomes.add(long_Chromosome.substring(i, i+30));
			
		}

		
		
		for(int i=0; i<crossedChromosomes.size()-1;i+=2){
		treeMap.put((float) ((-2*Math.pow(get_float(crossedChromosomes.get(i+1)), 3)+6*Math.pow(get_float(crossedChromosomes.get(i+1)), 2)+10)*Math.sin(Math.log(get_float(crossedChromosomes.get(i)))*Math.exp(get_float(crossedChromosomes.get(i+1)))))	, i)	; 
		}

		if(best_rezult[0]<treeMap.firstKey()){
		best_rezult[0]=treeMap.firstKey();
		best_rezult[1]=get_float(crossedChromosomes.get(treeMap.get(treeMap.firstKey())));
		best_rezult[2]=get_float(crossedChromosomes.get(treeMap.get(treeMap.firstKey()))+1);
		}
		
		
		// ннова попул€ц≥€ 
		int it=0;
		for(Map.Entry<Float, Integer> entry : treeMap.entrySet()){
			if(it>count)break;
			s_instance[it][0]=crossedChromosomes.get(entry.getValue());
			s_instance[it][1]=crossedChromosomes.get(entry.getValue()+1);
			it+=2;
			
		}
	
		treeMap.clear();
		crossedChromosomes.clear();

//		System.out.println("Iteration #"+ii);
//		System.out.println("Best rezult="+best_rezult[0] + " x1=" + bes t_rezult[1] + " x2=" + best_rezult[2]);
		 Float[] bestRezult=getBestRezult();
         best.put(bestRezult[0], bestRezult);
		}
		
		for(Map.Entry<Float, Float[]> entry : best.entrySet()){
			
//			System.out.println("Rezult= "+entry.getKey()+" x1= "+entry.getValue()[1]+" x2= "+ entry.getValue()[2]);
			System.out.printf(" Rezult= %.2f "+" x1= %.2f"+" x2= %.2f%n" ,entry.getKey(),entry.getValue()[1],entry.getValue()[2]);
		}
System.out.println("Best Rezult:");
		System.out.printf(" Rezult= %.2f "+" x1= %.2f"+" x2= %.2f%n" ,best.lastKey(),best.get(best.lastKey())[1],best.get(best.lastKey())[2]);

		timer_end=System.nanoTime();

		System.out.println("Time"+String.format("%,12d", timer_end-timet_start) + " ns") ;




	}

	public static float[][] get_FloatInstance (int n, int amount, float x1_min, float x1_max, float x2_min, float x2_max){
		float instance[][] = new float[n][amount];

		//generate instance
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < amount; j++) {
				if (j == 0)
					instance[i][j] = Float.parseFloat(df.format( rand.nextFloat() * (x1_max - x1_min) + x1_min )) ;
				if (j == 1)
					instance[i][j] = Float.parseFloat(df.format( rand.nextFloat() * (x2_max - x2_min) + x2_min )) ;
			}
		}

		return instance;
	}

	public static String[][] get_StringInstance (int n, int amount, float x1_min, float x1_max, float x2_min, float x2_max){
		String instance[][] = new String[n][amount];

		//generate instance
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < amount; j++) {
				if (j == 0)
					instance[i][j] = df.format( rand.nextFloat() * (x1_max - x1_min) + x1_min ) ;
				if (j == 1)
					instance[i][j] = df.format( rand.nextFloat() * (x2_max - x2_min) + x2_min ) ;
			}
		}

		return instance;
	}

	public static String[][] get_StringBinaryInstance (int n, int amount, float x1_min, float x1_max, float x2_min, float x2_max){
		String instance[][] = new String[n][amount];

		//generate instance
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < amount; j++) {
				if (j == 0)
					// WAT? генеруЇм значенн€ в д≥апазаон≥; df.format обр≥заЇм лишн≥ знаки(точн≥сть) отримуЇм String(float); Float.parseFloat отримуЇм Float(String); Float.floatToRawIntBits поб≥тово переводимось в int(float); Integer.toBinaryString отримуЇм дв≥йковий код String->Int(Float)
					instance[i][j] = Integer.toBinaryString(Float.floatToRawIntBits(Float.parseFloat(df.format( rand.nextFloat() * (x1_max - x1_min) + x1_min ))));
				if (j == 1)
					instance[i][j] = Integer.toBinaryString(Float.floatToRawIntBits(Float.parseFloat(df.format( rand.nextFloat() * (x2_max - x2_min) + x2_min ))));
			}
		}

		return instance;
	}

	public static String get_SimilarParent(String chromosome1, String[][] s_instance, int count ){
		String hromosome2 = null;
		int max_value=Integer.MIN_VALUE, index = 0;
		for(int j=0; j<count; j++ ){
			int i=0,sum=0; 
			for(char h: chromosome1.toCharArray()){
				if(h==s_instance[j][1].charAt(i)) sum++;
				i++; 
			}
			if(max_value<sum) {index=j; max_value =sum;}
		}
		hromosome2=s_instance[index][1];
		return hromosome2;
	}

	public static ArrayList<String> get_crossedChromosomes(ArrayList<String> crossedChromosomes,String chromosome1, String chromosome2){
		StringBuffer new_chromosome1 = new StringBuffer(chromosome1);
		StringBuffer new_chromosome2 = new StringBuffer(chromosome2);
		// схрещуванн€  хромосом. отримаэм 2 нов≥
		for(int i=1; i<chromosome1.length();i++){	
			if(i%2!=0)	new_chromosome1.setCharAt(i, chromosome1.charAt(i)); if (i%2==0) new_chromosome1.setCharAt(i, chromosome2.charAt(i-1));
			if(i%2==0)	new_chromosome2.setCharAt(i, chromosome2.charAt(i)); if (i%2!=0) new_chromosome2.setCharAt(i, chromosome1.charAt(i-1));	
		}

		crossedChromosomes.add(new_chromosome1.toString());
		crossedChromosomes.add(new_chromosome2.toString());
		return crossedChromosomes;
	}

	public static float get_float(String crosedChromosome){
		
		return Float.intBitsToFloat((Integer.parseInt(crosedChromosome,2)));
	}

	
	public static Float[] getBestRezult (){
		float x1 = 0,x2=0, rezult, bestRezult=0, param1, param2;
		
		for (int i = 0; i < n; i++) {	
					x1 = Float.parseFloat(df.format( rand.nextFloat() * (x1_max - x1_min) + x1_min )) ;
					x2 = Float.parseFloat(df.format( rand.nextFloat() * (x2_max - x2_min) + x2_min )) ;
					rezult = (float) ((-2*Math.pow(x2, 3)+6*Math.pow(x2, 2)+10)*Math.sin(Math.log(x1)*Math.exp(x2))); 
					if(rezult>bestRezult) {bestRezult=rezult; param1=x1; param2=x2;}
		}
		
//		System.out.println("#"+iterations+" rezult="+bestRezult+" x1= "+ x1+" x2="+x2);
		Float[] f ={bestRezult,x1,x2};
		iterations++;
		return f;
	}
	
	public static void bestHromosome(float[] bestRezult){
		float bestHtomosome[][] = new float[3][50];
	
	}
	

}
// €кий масив б≥льший? одновим≥рний чи двовим≥рний? €кий швидше працюЇ?
//String formatedDouble = String.format("%.2f", d);
//double d=-0.2; //17.875
// System.out.println(Double.doubleToLongBits(d));
// System.out.println(Double.doubleToRawLongBits(d));
// System.out.println(Long.toBinaryString(Double.doubleToLongBits(d)));
// System.out.println(Long.toBinaryString(Double.doubleToRawLongBits(d)));