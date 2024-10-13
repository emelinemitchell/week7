#include <iostream>
#include <string>
using namespace std;

bool isValidBase(char base){
    switch(base){ //takes users input
        case 'A':
            return true;
        case 'C':
            return true;
        case 'G':
            return true;
        case 'T':
            return true;
        default:
            return false; //if not A,C,T,G return false
    }
}

bool isValidStrand(string strand){
    //check if string only has A, C, G, T
    //use isValidBase
    if(strand.length() == 0){ //empty strand = false
        return false;
    }

    for(size_t  i = 0; i < strand.length(); i++){
        isValidBase(strand[i]); //check each index of length of string
        if(!isValidBase(strand[i])){ //checks if i != A, C, T, G
            return false;
        }
    }
    return true; //return true if strand is valid characters
}

double strandSimilarity(string strand1, string strand2){
    //check if string are equal length
    if(strand1.length() != strand2.length()){
        return 0; //if not equal return 0, end
    }
    //for loop to check if each index matchs
    //if index match, counter ++
    double match = 0; 
    double length = strand1.length(); //variable for length of string
    for(size_t i = 0; i < strand1.length(); i++){ //strand1 = strand2 so length is same
        if(strand1[i] == strand2[i]){
            match++; //if = match++
        }
    }
    //calculate similarity: similarVar / strand1.length()
    double similarity = match / length;
    return similarity;
}

int bestStrandMatch(string input_strand, string target_strand){
    size_t inputLength = input_strand.length();
    size_t targetLength = target_strand.length();

    if(targetLength > inputLength){ //returns -1 if target > input
        cout << "Best similarity score: 0.0" << endl; //score of 0
        return -1;
    }

    double similarity = 0; //similarity score
    double maxSimilarity = 0; //to find maximum
    double bestIndex = 0; //used to find best index

    for (size_t i = 0; i < input_strand.length(); i++){
        string check; 
        check = input_strand.substr(i, targetLength); //taking the string of input[i] to target length to find highest match

        //finds similarity of from i
        for(size_t j = 0; j < target_strand.length(); j++){
            similarity = strandSimilarity(check, target_strand); //checking similarity between target and section of input
            if(similarity == 0){
                similarity = 0; //if no simularities set similarity = 0
            }
        }
        if(similarity > maxSimilarity){ //finds max by testing against 0 then current max
            maxSimilarity = similarity;
            bestIndex = i; //best index is assosiated with max(best similarity)
        }  
    }

    cout << "Best similarity score: " << maxSimilarity << endl; 

    //returning index of best match
    if(maxSimilarity == 0){
        return 0;
    }else{
        return bestIndex; //if similarity != 0 then return the bestIndex where i has highest match
    }
}

//places 'x' where similarity strand isn't
string transformTargetWithX(string input, string target){
    size_t bestIndex = bestStrandMatch(input, target); //find similarity 
    if (bestIndex == 0 && target.length() > input.length()) { //If no valid match and target is longer
        return target; // Return original target
    }

    string result(input.length(), 'x'); // Create a new string filled with 'x'
    result.replace(bestIndex, target.length(), target); //place target strand with x's around it
    return result;
}

void identifyMutations(string input_strand, string target_strand) {
    size_t inputLength = input_strand.length();
    size_t targetLength = target_strand.length();

    // Test if matching and no mutations
    if (input_strand == target_strand){ //strings are same
        size_t bestIndex = bestStrandMatch(input_strand, target_strand);
        cout << "Best alignment index: " << bestIndex << endl;
        cout << "No mutations found." << endl;
        return; // If there are no mutations program ends
    }

    // Handle substitutions when strings are the same length
    if (inputLength == targetLength && input_strand != target_strand) {
        size_t bestIndex = bestStrandMatch(input_strand, target_strand);
        cout << "Best alignment index: " << bestIndex << endl;

        for (size_t i = 0; i < inputLength; i++) {
            if (input_strand[i] != target_strand[i]){ //if chars are != subsitution
                cout << "Substitution at position " << (i + 1) << ": " << input_strand[i] << " -> " << target_strand[i] << endl;
            }
        }
        return;
    }

    // Changes depending on if input or target is bigger
    if (inputLength > targetLength) {
        //Input larger - deletion/subsitution
        string modifyTarget = transformTargetWithX(input_strand, target_strand); //puts x before/after best similarity match position

        size_t bestIndex = 0;
        for (size_t i = 0; i < modifyTarget.length(); i++) {
            if (modifyTarget[i] != 'x'){
                bestIndex = i;
                break; // When it hits a character != 'x' break, when similarity strand starts
            }
        }
        cout << "Best alignment index: " << bestIndex << endl;

        //Checks for deletion and substitution (not insertion because input > target)
        for (size_t i = 0; i < inputLength; i++) {
            if (modifyTarget[i] == 'x'){
                cout << "Deletion at position " << (i + 1) << ": " << input_strand[i] << " is deleted in target strand" << endl;
            }
            if ((modifyTarget[i] != 'x') && (input_strand[i] != modifyTarget[i])){ //not 'x' but not == chars
                cout << "Substitution at position " << (i + 1) << ": " << input_strand[i] << " -> " << modifyTarget[i] << endl;
            }
        }
    } else {
        //Target larger - insertion/subsitution
        //same operations as if but reverses input/target
        string modifyInput = transformTargetWithX(target_strand, input_strand); 

        size_t bestIndex = 0;
        for (size_t i = 0; i < modifyInput.length(); i++) {
            if (modifyInput[i] != 'x') {
                bestIndex = i;
                break;
            }
        }
        cout << "Best alignment index: " << bestIndex << endl;

        //Checks for insertion and substitution (no deletion because target > input)
        for (size_t i = 0; i < targetLength; i++) {
            if (modifyInput[i] == 'x') {
                cout << "Insertion at position " << (i + 1) << ": " << target_strand[i] << " is inserted in target strand" << endl;
            }
            if ((modifyInput[i] != 'x') && (target_strand[i] != modifyInput[i])) {
                cout << "Substitution at position " << (i + 1) << ": " << modifyInput[i] << " -> " << target_strand[i] << endl;
            }
        }
    }
}

void transcribeDNAtoRNA(string strand){
    //change every 'T' to 'U'

    size_t strandLength = strand.length(); //variable for the length of strand

    //runs through each char and if 'T' changes to 'U'
    for(size_t i = 0; i < strandLength; i++){
        if(strand[i] == 'T'){
            strand[i] = 'U'; //replaces any char T with a U
        }
    }
    
    //print new RNA strand
    cout << "The transcribed DNA is: " << strand << endl;
}

void reverseComplement(string strand){
    size_t strandLength = strand.length(); //variable for the length of strand

    //reversing the string
    //for loop that reverses the string
    for(size_t i = 0; i < strandLength / 2; i++){ //only swap half the values, so it doesn't swap back
        char temp = strand[i]; //hold the character
        strand[i] = strand[strandLength - 1 - i]; //length of string - 1 - i to get the new position (swapped)
        strand[strandLength - 1 - i] = temp;
    }

    //change: A -> T, T -> A, C -> G, G -> C
    //take strand and loop through each character
        //for strand[i] change depending on character
    //switches the characters but doesn't put them in the reverse order
    for(size_t i = 0; i < strandLength; i++){
        //switch statement
        switch(strand[i]){
            case 'A':
                strand[i] = 'T';
                break; //so it doesn't switch back to A
            case 'T':
                strand[i] = 'A';
                break;
            case 'G':
                strand[i] = 'C';
                break;
            case 'C':
                strand[i] = 'G';
                break;
            //no default b/c we assume strand is all valid characters 
        }
    }
    cout << "The reverse complement is: " << strand << endl; //just print the reversed string
} 

void getCodingFrames(string strand){
    size_t strandLength = strand.length();
    //for loop to go all the way through one at a time
    //find the first atg then stop

    size_t frameCounter = 0; 
    size_t starting_index = 0;
    while (starting_index < strandLength){ //index goes up by 1, loop continues for the length of string
        bool start = false;

        //find start codon
        for (size_t i = starting_index; i < strandLength - 1; i++) {
            if(strand.substr(i, 3) == "ATG"){ //if ATG make start true and index + 1
                starting_index = i;
                start = true;
                break; //exit for loop when find ATG
            }
        }

        //find stop codon
        bool stop_codon_found = false;
        if (start){ 
        for(size_t i = starting_index; i < strandLength - 2; i += 3){ 
            string current_codon = strand.substr(i, 3);

            if(current_codon == "TAA" || current_codon == "TAG" || current_codon == "TGA"){
                //stoping codon
                cout << "The extracted reading frames are: " << endl;
                cout << strand.substr(starting_index, (i + 3) - starting_index) << endl;
                frameCounter++;
                stop_codon_found = true;
                start = false;
                starting_index = i;
                break;
            }
        }
        }
        if (!stop_codon_found) {
            starting_index += 1;
        }
        //still need to handle if you dont have reading frames
    }
    if(frameCounter == 0){ //if not ATG or ending codon there are no frames 
        cout << "The extracted reading frames are: " << endl;
        cout << "No reading frames found." << endl;
    }
    return;
}

int main(){
    size_t selector = 0;

    while(selector != 7){
        cout << "--- DNA Analysis Menu ---" << endl;
        cout << "1. Calculate the similarity between two sequences of the same length" << endl;
        cout << "2. Calculate the best similarity between two sequences of either equal or unequal length" << endl;
        cout << "3. Identify mutations" << endl;
        cout << "4. Transcribe DNA to RNA" << endl;
        cout << "5. Find the reverse complement of a DNA sequence" << endl;
        cout << "6. Extract coding frames" << endl;
        cout << "7. Exit" << endl;
        cout << "Please enter your choice (1 - 7):" << endl;

        string strand1 = ""; //initalizes strand1 and strand2 to empty to be called
        string strand2 = "";
        double similarity = 0.0; //initalize outside switch for compiling errors, simplicity 

        cin >> selector;

        switch(selector){
            case 1: //calculate similarity of sequences of = length
                cout << "Enter the first DNA sequence: " << endl;
                cin >> strand1;
                while(!isValidStrand(strand1)){ //if not a valid strand ask for new input
                    cout << "Invalid input. Please enter a valid sequence." << endl;
                    cout << "Enter the first DNA sequence: " << endl;
                    cin >> strand1;
                }

                cout << "Enter the second DNA sequence: " << endl;
                cin >> strand2;
                while(!isValidStrand(strand2)){
                    cout << "Invalid input. Please enter a valid sequence." << endl;
                    cout << "Enter the second DNA sequence: " << endl;
                    cin >> strand2;
                }


                //if strands aren't valid or != exit case
                if(!isValidStrand(strand1) && !isValidStrand(strand2)){
                    cout << "Error: Input strands must be of the same length." << endl;
                    break;
                }
                if(strand1.length() != strand2.length()){
                    cout << "Error: Input strands must be of the same length." << endl;
                    break;
                }


                similarity = strandSimilarity(strand1, strand2); //calculates similarity and prints
                cout << "Similarity score: " << similarity << endl;
                break;

            case 2: //calculate similarity of sequences of = or != length
                cout << "Enter the first DNA sequence: " << endl;
                cin >> strand1;
                while(!isValidStrand(strand1)){
                    cout << "Invalid input. Please enter a valid sequence." << endl;
                    cout << "Enter the first DNA sequence: " << endl;
                    cin >> strand1;
                }

                cout << "Enter the second DNA sequence: " << endl;
                cin >> strand2;
                while(!isValidStrand(strand2)){
                    cout << "Invalid input. Please enter a valid sequence." << endl;
                    cout << "Enter the second DNA sequence: " << endl;
                    cin >> strand2;
                }

                bestStrandMatch(strand1, strand2);
                break;

            case 3: //Identify Mutations
                cout << "Enter the first DNA sequence: " << endl;
                cin >> strand1;
                while(!isValidStrand(strand1)){
                    cout << "Invalid input. Please enter a valid sequence." << endl;
                    cout << "Enter the first DNA sequence: " << endl;
                    cin >> strand1;
                }

                cout << "Enter the second DNA sequence: " << endl;
                cin >> strand2;
                while(!isValidStrand(strand2)){
                    cout << "Invalid input. Please enter a valid sequence." << endl;
                    cout << "Enter the first DNA sequence: " << endl;
                    cin >> strand2;
                }

                //funtion will print out mutations, similarity, and alignment
                identifyMutations(strand1, strand2);
                break;

            case 4: //Transcribe DNA to RNA
                cout << "Enter the DNA sequence to be transcribed: " << endl;
                cin >> strand1;
                while(!isValidStrand(strand1)){
                    cout << "Invalid input. Please enter a valid sequence." << endl;
                    cout << "Enter the DNA sequence to be transcribed: " << endl;
                    cin >> strand1;
                }

                //print the transcribed DNA in function
                transcribeDNAtoRNA(strand1);
                break;

            case 5: //Reverse complement of DNA sequence 
                cout << "Enter the DNA sequence: " << endl;
                cin >> strand1;
                while(!isValidStrand(strand1)){
                    cout << "Invalid input. Please enter a valid sequence." << endl;
                    cout << "Enter the DNA sequence: " << endl;
                    cin >> strand1;
                }

                //function will print out reversed strand
                reverseComplement(strand1);
                break;

            case 6: //extract reading frames
                cout << "Enter the DNA sequence: " << endl;
                cin >> strand1;
                while(!isValidStrand(strand1)){
                    cout << "Invalid input. Please enter a valid sequence." << endl;
                    cout << "Enter the DNA sequence: " << endl;
                    cin >> strand1;
                }

                //function will print out frames
                getCodingFrames(strand1);
                break;
                
            case 7: //exit 
                cout << "Exiting program." << endl;
                break;

            default:
                cout << "Invalid input. Please select a valid option." << endl; //loop starts over 
                break;
        }
    }
    return 0;
}
