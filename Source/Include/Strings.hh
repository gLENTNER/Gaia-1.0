// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved
// See LICENSE file (GPLv3)
// Include/Strings.hh
//
// Helper functions for manipulating strings

#ifndef _STRINGS_HH_
#define _STRINGS_HH_

namespace GAIA {
	
// remove all characters after `delim`
void Clip(std::string &input_string, const std::string &delim){
	
	std::size_t pos = input_string.find(delim);
	
	if ( pos != std::string::npos )
		input_string.replace(pos, input_string.length(), "");
}

// split a string into `words`, grouping quoted words
std::vector<std::string> Split(const std::string &input){

	std::vector<std::string> sentence;
	std::string word;
	bool quoted = false;

	for (int i = 0; i < input.length(); i++){
	
		if (input[i] == ' '){
			
			if (quoted) word += " ";
			else {
			
				if ( !word.empty() ) 
					sentence.push_back(word);
				
				word = "";
			}
		
		}

		else if ( input[i] == '\"' ) quoted = !quoted;
		else word += input[i];
	}
	
	// get last word
	if ( !word.empty() )
		sentence.push_back(word);

	return sentence;
}

// replace all instances of `search_str` in `input_str` with `replace_str`
void ReplaceAll(const std::string &search_str, 
	const std::string &replace_str, std::string& input_str ){
	
	std::size_t pos = 0; 
	
	while ( ( pos = input_str.find(search_str, pos) ) != std::string::npos ){
	
		input_str.replace( pos, search_str.length( ), replace_str );
		pos += replace_str.length( );
	}
}

} // namespace GAIA

#endif