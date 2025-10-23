#include <Rcpp.h>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <vector>
#include <sstream>
//[[Rcpp::depends(BH)]]
//#include <boost/unordered/unordered_map.hpp>

using namespace Rcpp;
using namespace std;

// Helper function to split a string based on a delimiter
vector<string> split_string(const string& str, const string& delimiter) {
    vector<string> tokens;
    size_t start = 0, end = 0;
    while ((end = str.find(delimiter, start)) != string::npos) {
        tokens.push_back(str.substr(start, end - start));
        start = end + delimiter.length();
    }
    tokens.push_back(str.substr(start));
    return tokens;
}

// Helper function to join a vector of strings with a separator
string join_strings(const vector<string>& strings, const string& separator) {
    stringstream result;
    for (size_t i = 0; i < strings.size(); ++i) {
        if (i > 0) result << separator;
        result << strings[i];
    }
    return result.str();
}

// Function to escape JSON special characters
string escape_json_string(const string& str) {
    stringstream escaped;
    for (char c : str) {
        switch (c) {
            case '"': escaped << "\\\""; break;
            case '\\': escaped << "\\\\"; break;
            case '/': escaped << "\\/"; break;
            case '\b': escaped << "\\b"; break;
            case '\f': escaped << "\\f"; break;
            case '\n': escaped << "\\n"; break;
            case '\r': escaped << "\\r"; break;
            case '\t': escaped << "\\t"; break;
            default: escaped << c;
        }
    }
    return escaped.str();
}

// [[Rcpp::export]]
std::string paste_action_tokens_json_cpp(CharacterVector parameters, 
                                         CharacterVector values, 
                                         CharacterVector fullstrings, 
                                         DataFrame labactions, 
                                         std::string split = "\\^") {
    
    // Convert R input to C++ data structures
    int n = fullstrings.size();
    vector<vector<string>> tokenlist(n);
    unordered_set<string> uniquetokens;
    
    // Split fullstrings and build the set of unique tokens
    for (int i = 0; i < n; ++i) {
        tokenlist[i] = split_string(as<string>(fullstrings[i]), split);
        for (const auto& token : tokenlist[i]) {
            uniquetokens.insert(token);
        }
    }
    
    // Extract labactions columns
    CharacterVector action = labactions["action"];
    CharacterVector statement = labactions["Statement"];
    int labactions_size = statement.size();
    
    // Map from statement to action for easy lookup
    unordered_map<string, string> statement_action_map;
	//Heh, use boost::unordered_flat_map for fun :)
	//boost::unordered_flat_map<string, string> statement_action_map;
    for (int i = 0; i < labactions_size; ++i) {
        if (uniquetokens.count(as<string>(statement[i])) > 0) {
            statement_action_map[as<string>(statement[i])] = as<string>(action[i]);
        }
    }
    
    // Build the JSON result
    stringstream json_output;
    json_output << "{ \"actions\": [";
    
    bool first_action = true;
    for (const auto& [stmt, act] : statement_action_map) {
        if (!first_action) json_output << ", ";
        first_action = false;
        
        // Build affected and deleted parameters
        vector<string> affected_params;
        vector<string> deleted_params;
        
        for (int i = 0; i < n; ++i) {
            if (find(tokenlist[i].begin(), tokenlist[i].end(), stmt) != tokenlist[i].end()) {
                if (act != "Delete") {
                    affected_params.push_back(escape_json_string(as<string>(parameters[i])));
                } else {
                    string del_entry = "{\"" + 
						escape_json_string(as<string>(parameters[i])) + 
						"\" : " + 
						escape_json_string(as<string>(values[i])) + 
						"}";
                    deleted_params.push_back(del_entry);
                }
            }
        }
        
        // Create JSON structure for the current action
        json_output << "{";
        json_output << "\"st\": \"" << escape_json_string(stmt) << "\", ";
        json_output << "\"act\": \"" << escape_json_string(act) << "\", ";
        json_output << "\"cols\": {";
        if (act == "Delete") {
            json_output << "\"deleted\": [" << join_strings(deleted_params, ", ") << "]";
        } else {
            json_output << "\"affected\": [\"" << join_strings(affected_params, "\", \"") << "\"]";
        }
        json_output << "} }";
    }
    
    json_output << "] }";
    
    return json_output.str();
}