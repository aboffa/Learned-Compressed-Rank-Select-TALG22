// This file is part of la_vector <https://github.com/aboffa/Learned-Compressed-Rank-Select-TALG22>.
// Copyright (c) 2022 Antonio Boffa.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef IMPROVING_LA_VECTOR_MY_UTILS_H
#define IMPROVING_LA_VECTOR_MY_UTILS_H

class error_occurrence_counter{
private:
    std::map<int,int> occurrences;
public:
    void add(int error){
        auto it = occurrences.find(error);
        if(it != occurrences.end())
            it->second++;
        else
            occurrences[error] = 1;

    }

    std::map<int,int> get(){
        std::map<int,int> to_return(occurrences);
        clear();
        return to_return;
    }

    void clear(){
        occurrences.clear();
    }
};
// Global object used to count the number of occurrences the corrections
error_occurrence_counter eoc;

class u_corrections_vec{
private:
    std::vector<uint16_t> u_cvec;
public:
    void add(int error){
        u_cvec.push_back(error);
    }
    void clear(){
        u_cvec.clear();
    }

    std::vector<uint16_t> get(){
        std::vector<uint16_t> to_return(u_cvec);
        clear();
        return to_return;
    }
};

// Global object used to collect the corrections
u_corrections_vec u_cvec;

#endif //IMPROVING_LA_VECTOR_MY_UTILS_H
