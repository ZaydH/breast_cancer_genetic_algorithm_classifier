package cs123A;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class BreastCancerDataSet {

	private static int trainingDataSetSize = 200;
	public final static int MAXIMUM_TRAINING_DATA_SET_SIZE = 682;
	private List<Patient> setOfPatients;
	
	private static boolean allowShuffling = true;
	
	/**
	 * Constructor to build datasets of breast cancer dataset.
	 */
	public BreastCancerDataSet(){
		//----- Create the patient array
		setOfPatients = new ArrayList<Patient>();
	}
	
	/**
	 * Additional Private Constructor that takes an ArrayList of patients.
	 * 
	 * @param patients ArrayList of patients that will serve as the patients in the data set.
	 */
	private BreastCancerDataSet(List<Patient> patients){
		setOfPatients = patients;
	}
	
	
	/**
	 * Accessor for the Training Data Set Size
	 * 
	 * @return Size of the Training Data Set
	 */
	public static int getTrainingDataSetSize(){
		return trainingDataSetSize;
	}
	

	/**
	 * Mutator for the Training Data Set Size.
	 */
	public static void setTrainingDataSetSize(int newTrainingDataSetSize){
		trainingDataSetSize = newTrainingDataSetSize;
	}
	
	
	/**
	 * Adds a new patient to the training set.
	 * 
	 * @param features A comma separated list of the features that is parsed by the patient.
	 */
	public void addPatient(String features){
		Patient newPatient = new Patient(features);
		setOfPatients.add(newPatient);
	}
	

	/**
	 * Adds a new patient to the training set.
	 * 
	 * @param newPatient A Patient object,
	 */
	public void addPatient(Patient newPatient){
		setOfPatients.add(newPatient);
	}
	
	
	/**
	 * Accessor for number of patients in the dataset.
	 * 
	 * @return Number of patients in the data set.
	 */
	public int getDataSetSize(){
		return setOfPatients.size();
	}
	
	
	/**
	 * Splits a breast cancer into two.  The implicit data set is reduced in size by
	 * "numbElements" and those elements are placed into the returned BreastCancerDataSet
	 * object.
	 * 
	 * @param numbElements Number of elements to remove from the current data
	 * @return New BreastCancerDataSet of size numbElements
	 */
	public BreastCancerDataSet removeRandomSubset(int numbElements){
		//---- Shuffle the ArrayList.
		if(allowShuffling)
			Collections.shuffle(setOfPatients);
		
		List<Patient> removedSublist = new ArrayList<Patient>(setOfPatients.subList(0, numbElements));
		
		//----- Update this objects set of patients.
		setOfPatients = new ArrayList<Patient>(setOfPatients.subList(numbElements, setOfPatients.size()));
		
		//---- Create a new breast cancer data set.
		return new BreastCancerDataSet(removedSublist);
	}
	
	
	/**
	 * This is used to disable shuffling in the population for repeatibility analysis
	 */
	public static void disableDataSetRandomShuffling(){
		allowShuffling = false;
	}
	
	
	/**
	 * Static method to merge two BreastCancerDataSets.
	 * 
	 * @param dataSet1 First BreastCancerDataSet object.
	 * @param dataSet2 Second BreastCancerDataSet object.
	 * 
	 * @return Merged  BreastCancerDataSet object.
	 */
	public static BreastCancerDataSet mergeDataSets(BreastCancerDataSet dataSet1, 
													BreastCancerDataSet dataSet2){
		
		//--- Combine the two lists of the datasets.
		List<Patient> mergedList = new ArrayList<Patient>();
		
		//---- Append the two lists.
		mergedList.addAll(dataSet1.setOfPatients);
		mergedList.addAll(dataSet2.setOfPatients);
		
		//--- Build and return the new data set.
		return new BreastCancerDataSet(mergedList);
		
	}
	
	
	
	/**
	 * Determines the score of a given chromosome and a population.
	 * 
	 * @param chromosome Chromosome whose population score will be calculated.
	 * @return Score for the chromosome which is the number of correct identifications versus incorrect identifications.
	 */
	public int getChromosomeScoreForPopulation(GAChromosome chromosome){
		return getChromosomeScoreForPopulation(chromosome, 1);
	}
	
	/**
	 * Determines the score of a given chromosome and a population.
	 * 
	 * @param chromosome Chromosome whose population score will be calculated.
	 * @param malignancyBiasFactor A bias factor to skew the results to favor correct scoring of malignant tumors.
	 * @return Score for the chromosome which is the number of correct identifications versus incorrect identifications.
	 */
	public int getChromosomeScoreForPopulation(GAChromosome chromosome, int malignancyBiasFactor){
		return (int)Math.round(getChromosomeScoreAndSeparationForPopulation(chromosome, 1)[0]);
	}
	

	public double[] getChromosomeScoreAndSeparationForPopulation(GAChromosome chromosome, int malignancyBiasFactor){
		
		int chromosomeScore = 0;
		int index;
		double patientScore;
		Patient patient;
		double separation = 0;
		
		int[] chromosomeGainVector = chromosome.getGainVector();
		int chromosomeOffset = chromosome.getOffset();
		
		//---- Iterate through the 
		for(index = 0; index < setOfPatients.size(); index++){
			//---- Get the patient
			patient = setOfPatients.get(index);
			//---- Determine the score for that patient
			patientScore = patient.calculateLinearFunction(chromosomeGainVector, chromosomeOffset);
			//---- Update the separation.
			separation += patientScore;
			//---- If the patient is correctly categorized, give it a positive score.
			if(patientScore > 0)
				chromosomeScore++;
				if(patient.isMalignant())
					chromosomeScore += (malignancyBiasFactor-1); //--- Subtract one since already incremented
			else
				if(patient.isMalignant())
					chromosomeScore -= (malignancyBiasFactor-1); //--- Subtract one since already decremented
//			//---- Patient is miscategorized so give it a negative score.
//			else
//				chromosomeScore--;	
			
		}
		
		//---- Give the chromosome score.
		return new double[] {chromosomeScore, separation};
		
	}
	
	/**
	 * Calculates the accuracy for malignant tumors properly categorized.
	 * 
	 * @param chromosome	Genetic Algorithm Chromosome (Linear function)
	 * @return				Percentage of malignant tumors correctly categorized.
	 */
	public double getMaligancyAccuracyForPopulation(GAChromosome chromosome){
		
		int index;
		int malignancyCorrect = 0;
		int numberOfMalignantPatients = 0;
		long patientScore;
		Patient patient;
		
		
		int[] chromosomeGainVector = chromosome.getGainVector();
		int chromosomeOffset = chromosome.getOffset();
		
		//---- Iterate through the 
		for(index = 0; index < setOfPatients.size(); index++){
			//---- Get the patient
			patient = setOfPatients.get(index);
			if(patient.isMalignant()){
				numberOfMalignantPatients++;
				
				//---- Determine the score for that patient
				patientScore = patient.calculateLinearFunction(chromosomeGainVector, chromosomeOffset);
				//---- If the patient is correctly categorized, give it a positive score.
				if(patientScore > 0)
					malignancyCorrect++;
			}	
		}
		
		//---- Give the chromosome score.
		return malignancyCorrect * 100.0 / numberOfMalignantPatients ;
		
	}
	
	
	/**
	 * Determines the percentage of correct assignments for this data set and chromosome pairing.
	 * 
	 * @param chromosome GAChromosome Object
	 * @return Percent of correct classifications 
	 */
	public double getPercentCorrect(GAChromosome chromosome){
		int numbCorrect = this.getChromosomeScoreForPopulation(chromosome);
		return numbCorrect * 100.0 / this.getDataSetSize();
	}
	
}
