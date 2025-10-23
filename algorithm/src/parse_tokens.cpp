#include <Rcpp.h>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <vector>
#include <sstream>

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

// [[Rcpp::export]]
std::string paste_action_tokens_cpp(CharacterVector parameters, 
                                    CharacterVector values, 
                                    CharacterVector fullstrings, 
                                    DataFrame labactions, 
                                    std::string split = "^") {
    
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
    for (int i = 0; i < labactions_size; ++i) {
        if (uniquetokens.count(as<string>(statement[i])) > 0) {
            statement_action_map[as<string>(statement[i])] = as<string>(action[i]);
        }
    }
    
    // Build the results for affected and deleted columns
    vector<string> cols_affected, cols_deleted;
    
    for (const auto& [stmt, act] : statement_action_map) {
        vector<string> affected_params;
        vector<string> deleted_params;
        
        // Check each tokenlist for the presence of the statement
        for (int i = 0; i < n; ++i) {
            if (find(tokenlist[i].begin(), tokenlist[i].end(), stmt) != tokenlist[i].end()) {
                if (act != "Delete") {
                    affected_params.push_back(as<string>(parameters[i]));
                } else {
                    string del_entry = as<string>(parameters[i]) + " (" + as<string>(values[i]) + ")";
                    deleted_params.push_back(del_entry);
                }
            }
        }
        
        // Join results
        if (act != "Delete" && !affected_params.empty()) {
            cols_affected.push_back(join_strings(affected_params, ","));
        }
        if (act == "Delete" && !deleted_params.empty()) {
            cols_deleted.push_back(join_strings(deleted_params, ", "));
        }
    }
    
    // Combine affected and deleted columns
    vector<string> allcols;
    allcols.insert(allcols.end(), cols_deleted.begin(), cols_deleted.end());
    allcols.insert(allcols.end(), cols_affected.begin(), cols_affected.end());
    
    // Construct final result
    vector<string> final_result;
    for (const auto& [stmt, act] : statement_action_map) {
        string result = stmt + " (" + act + ") : [ " + (act == "Delete" ? cols_deleted[0] : cols_affected[0]) + " ]";
        final_result.push_back(result);
    }
    
    return join_strings(final_result, ", ");
}

// [[Rcpp::export]]
NumericVector mark_deleted_values_cpp(NumericVector values, 
                                      CharacterVector fullstrings, 
                                      DataFrame labactions, 
                                      std::string split = "^") {
    
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
    
    // Create a set of tokens that should be deleted
    unordered_set<string> delete_tokens;
    for (int i = 0; i < labactions_size; ++i) {
        if (as<string>(action[i]) == "Delete" && uniquetokens.count(as<string>(statement[i])) > 0) {
            delete_tokens.insert(as<string>(statement[i]));
        }
    }
    
    // Create a copy of values to modify
    NumericVector modified_values = clone(values);
    
    // Check each tokenlist for tokens that are in the delete set
    for (int i = 0; i < n; ++i) {
        for (const auto& token : tokenlist[i]) {
            if (delete_tokens.count(token) > 0) {
                // If a delete token is present, set the corresponding value to NA
                modified_values[i] = NA_REAL;
                break; // No need to check further tokens for this row
            }
        }
    }
    
    return modified_values;
}