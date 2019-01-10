//~ de Bruijn graph TRIMing tool
//~ Copyright (C) 2017  Limasset Antoine and Marchet Camille

//~ This program is free software: you can redistribute it and/or modify
//~ it under the terms of the GNU Affero General Public License as
//~ published by the Free Software Foundation, either version 3 of the
//~ License, or (at your option) any later version.

//~ This program is distributed in the hope that it will be useful,
//~ but WITHOUT ANY WARRANTY; without even the implied warranty of
//~ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//~ GNU Affero General Public License for more details.

//~ You should have received a copy of the GNU Affero General Public License
//~ along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <fstream>
#include <cstring>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <bitset>
#include <chrono>
#include <unistd.h>
#include <stdio.h>





using namespace std;
using namespace chrono;




string intToString(uint64_t n){
	if(n<1000){
		return to_string(n);
	}
	string end(to_string(n%1000));
	if(end.size()==3){
		return intToString(n/1000)+","+end;
	}
	if(end.size()==2){
		return intToString(n/1000)+",0"+end;
	}
	return intToString(n/1000)+",00"+end;
}



char revCompChar(char c) {
	switch (c) {
		case 'A': return 'T';
		case 'C': return 'G';
		case 'G': return 'C';
	}
	return 'A';
}



string revComp(const string& s){
	string rc(s.size(),0);
	for (int i((int)s.length() - 1); i >= 0; i--){
		rc[s.size()-1-i] = revCompChar(s[i]);
	}
	return rc;
}


string get_canon(const string& s){
	return min(s,revComp(s));
}



uint64_t str2num(const string& str){
	uint64_t res(0);
	for(uint i(0);i<str.size();i++){
		res<<=2;
		switch (str[i]){
			case 'A':res+=0;break;
			case 'C':res+=1;break;
			case 'G':res+=2;break;
			default:res+=3;break;
		}
	}
	return res;
}



vector<bool> str2bool(const string& str){
	vector<bool> result;
	for(uint i(0);i<str.size();i++){
		switch (str[i]){
			case 'A':result.push_back(false);result.push_back(false);break;
			case 'C':result.push_back(false);result.push_back(true);break;
			case 'G':result.push_back(true);result.push_back(false);break;
			default:result.push_back(true);result.push_back(true);break;
		}
	}
	return result;
}



vector<bool> int2bool(uint n){
	vector<bool> result;
	while(n>0){
		if(n%2==0){
			result.push_back(false);
		}else{
			result.push_back(true);
		}
		n>>=1;
	}
	if(result.size()%2==0){
		result.push_back(false);
	}
	return result;
}



uint bool2int(const vector<bool>& v){
	uint result(0);
	uint acc(1);
	for(uint i(0);i<v.size();++i){
		if(v[i]){
			result+=acc;
		}
		acc<<=1;
	}
	return result;
}



string bool2str(const vector<bool>& v){
	string result;
	for(uint i(0);i<v.size();i+=2){
		if(v[i]){
			if(v[i+1]){
				result.push_back('T');
			}else{
				result.push_back('G');
			}
		}else{
			if(v[i+1]){
				result.push_back('C');
			}else{
				result.push_back('A');
			}
		}
	}
	return result;
}



uint32_t xs(uint32_t y){
	y^=(y<<13); y^=(y>>17); return (y^=(y<<15));
}



char randNucle(char c){
	switch (rand()%4){
		case 0:
			if(c!='A'){
				return 'A';
			}
		case 1:
			if(c!='C'){
				return 'C';
			}
		case 2:
			if(c!='G'){
				return 'G';
			}
		case 3:
			if(c!='T'){
				return 'T';
			}
	}
	return randNucle(c);
}



string compactionNoRecur(const string& seq1,const string& seq2, uint k){
	//~ cout<<seq1<<" "<<seq2<<endl;
	uint s1(seq1.size()),s2(seq2.size());
	if(s1==0 or s2==0){return "";}
	string rc2(revComp(seq2)),end1(seq1.substr(s1-k,k)), beg2(seq2.substr(0,k));
	if(end1==beg2){return seq1+(seq2.substr(k));}
	string begrc2(rc2.substr(0,k));
	if(end1==begrc2){return seq1+(rc2.substr(k));}
	//~ return compaction(revComp(seq1),seq2,k);
	cout<<"fail"<<endl;
	return "";
}


string compaction(const string& seq1,const string& seq2, uint k){
	//~ cout<<seq1<<" "<<seq2<<endl;
	uint s1(seq1.size()),s2(seq2.size());
	if(s1==0 or s2==0){return "";}
	string rc2(revComp(seq2)),end1(seq1.substr(s1-k,k)), beg2(seq2.substr(0,k));
	if(end1==beg2){return seq1+(seq2.substr(k));}
	string begrc2(rc2.substr(0,k));
	if(end1==begrc2){return seq1+(rc2.substr(k));}
	return compactionNoRecur(revComp(seq1),seq2,k);
}



uint getPosition(const vector<vector<bool>>& unitigs,uint n){
	//~ cout<<"getPOsition"<<n<<endl;
	auto v(unitigs[n]);
	if(v.size()%2==1){
		uint newPos(bool2int(v));
		return getPosition(unitigs,newPos);
	}else{
		return n;
	}
	return 0;
}



void str2bin(const string& str, char* cstr){
	uint j(0),mult(1);
	for(uint i(0); i<str.size(); ++i){
		switch (str[i]){
			//~ case 'A':cstr[j]+=0;break;
			case 'C':cstr[j]+=1*mult;break;
			case 'G':cstr[j]+=2*mult;break;
			case 'T':cstr[j]+=3*mult;break;
		}
		mult<<=2;
		if(mult>64){
			++j;
			mult=1;
		}
	}
}



vector<pair<string,uint32_t>> uniqueOnly(vector<pair<string,uint32_t>>& v){
	bool unique(true);
	vector<pair<string,uint32_t>> result;
	for(uint i(0);i+1<v.size();++i){
		if(v[i].first==v[i+1].first){
			unique=false;
		}else{
			if(unique){
				result.push_back(v[i]);
			}
			unique=true;
		}
	}
	if(unique and v.size()>0){
		result.push_back(v[v.size()-1]);
	}
	return result;
}



double parseCoverage(const string& str){
	size_t pos(str.find("km:f:"));
	if(pos==string::npos){
		pos=(str.find("KM:f:"));
	}
	if(pos==string::npos){
		return 1;
	}
	uint i(1);
	while(str[i+pos+5]!=' '){
		++i;
	}
	return stof(str.substr(pos+5,i));
}



void usage(){
	cout<<"Usage:"<<endl;
		cout<<"-u [Unitig file]\n"
		<<"-k [Kmer size]\n"
		<<"-t [Tipping length (none)]\n"
		<<"-T [Cleaning Step (1)]\n"
		<<"-c [Core used (1)]\n"
		<<"-h [Hash size, use 2^h files (8 for 256 files)]\n"
		<<"-f [Unitig min coverage (none, 0 for auto)]\n"
		<<"-m [Unitig max coverage (none)]\n"
		<<"-a [Edge filtering ratio (none)]\n"
		<<"-l low coverage\n"
		<<"-l high coverage\n"
		<<"-o [Output file (out_tipped)]\n"
		<<endl;
}



struct args_btt{
	uint nbFiles;
	uint coreUsed;
	uint kmerSize;
	uint unitigThreshold;
	uint unitigThreshold_MAX;
	uint tipingSize;
	uint ratioCoverage;
	uint low_coverage;
	uint high_coverage;
};


bool parallel_unitigs(const vector<bool>& u1,const vector<bool> u2, uint k){
	if(u1.size()!=u2.size()){
		return false;
	}
	//~ cout<<"go"<<endl;
	vector<bool> begin1 (u1.begin(), u1.begin()+2*k);
	vector<bool> begin2(u2.begin(), u2.begin()+2*k);
	string canon_begin1(get_canon(bool2str(begin1)));
	string canon_begin2(get_canon(bool2str(begin2)));
	vector<bool> end1( u1.end()-2*k, u1.end());
	vector<bool> end2( u2.end()-2*k, u2.end());
	string canon_end1(get_canon(bool2str(end1)));
	string canon_end2(get_canon(bool2str(end2)));
	//~ if(not ((canon_begin1==canon_begin2 or canon_begin1== canon_end2) and (canon_end1==canon_begin2 or canon_end1==canon_end2))){
		//~ cout<<bool2str(u1)<<endl;
		//~ cout<<bool2str(u2)<<endl;
		//~ cout<<canon_end1<<" "<<canon_end2<<endl;
		//~ cout<<canon_begin1<<" "<<canon_begin2<<endl;

			//~ cout<<"fAIL"<<endl;

		//~ cin.get();
	//~ }
	return ((canon_begin1==canon_begin2 or canon_begin1== canon_end2) and (canon_end1==canon_begin2 or canon_end1==canon_end2));
}

bool less_than_n_missmatch(const string& str1, const string& str2, uint max_miss){
	uint min_size(min(str1.size(),str2.size()));
	uint miss(0);
	for(uint i(0);i<min_size;++i){
		if(str1[i]!=str2[i]){
			if(++miss>max_miss){
				return false;
			}
		}
	}
	return true;
}


bool low_sistance_tip(const vector<bool>& Btip,vector<bool>& Bref, const string& overlap){
	//~ cout<<"go"<<endl;
	uint k=overlap.size();
	string tip(bool2str(Btip));
	if(get_canon(tip.substr(0,k))!=overlap){
		tip=revComp(tip);
	}
	string ref(bool2str(Bref));
	if(get_canon(ref.substr(0,k))!=overlap){
		ref=revComp(ref);
	}
	//~ cout<<overlap<<endl;
	//~ cout<<tip<<"\n"<<ref<<endl;cin.get();
	//~ if(less_than_n_missmatch(tip,ref,2)){
		//~ cout<<tip<<endl;
		//~ cout<<ref<<endl;cin.get();
	//~ }
	return less_than_n_missmatch(tip,ref,1);
}



int cleaning(string outFile, string inputUnitig,args_btt arg){
	auto start=system_clock::now();
	uint64_t tiping(0),compactions(0),unitigFiltered(0),advanced_tipping(0),island(0),bulles(0);
	ofstream out(outFile);

	ifstream inUnitigs(inputUnitig);
	vector<fstream> beginFiles(arg.nbFiles),endFiles(arg.nbFiles);
	for(uint i(0); i< arg.nbFiles; ++i){
		beginFiles[i].open(".begin"+to_string(i),fstream::out|fstream::in|fstream::binary|fstream::trunc);
		endFiles[i].open(".end"+to_string(i),fstream::out|fstream::in|fstream::binary|fstream::trunc);
	}
	//THIS MEMORY USAGE CAN BE REMOVED
	vector<vector<bool>> unitigs;
	vector<uint> coverages;



	cout<<"\tParitioning"<<endl;
	//PARTITIONING
	string begin,end,unitig,useless,beginRc,endRc;
	uint64_t hashBegin, hashBeginRc, hashEnd, hashEndRc;
	uint32_t unitigIndice(0);
	vector<uint> abundance_unitigs(100);
	while(not inUnitigs.eof()){
		getline(inUnitigs,useless);
		unitig="";
		getline(inUnitigs,unitig);
		if(unitig.size()<arg.kmerSize){
			continue;
		}
		uint coverage=parseCoverage(useless);
		if(((arg.unitigThreshold>1 and coverage<arg.unitigThreshold)) or (coverage>arg.unitigThreshold_MAX and arg.unitigThreshold_MAX>0)){
			unitigFiltered++;
			continue;
		}
		unitigs.push_back(str2bool(unitig));
		//~ cout<<coverage<<endl;
		coverages.push_back(coverage);
		if(coverage<100){
			abundance_unitigs[coverage]++;
		}else{
			abundance_unitigs[0]++;
		}


		begin=unitig.substr(0,arg.kmerSize);
		beginRc=revComp(begin);
		hashBegin=str2num(begin);
		hashBeginRc=str2num(beginRc);
		if(hashBegin<=hashBeginRc){
			beginFiles[xs(hashBegin)%arg.nbFiles]<<begin;
			beginFiles[xs(hashBegin)%arg.nbFiles].write(reinterpret_cast<const char*>(&unitigIndice), 4);
		}else{
			endFiles[xs(hashBeginRc)%arg.nbFiles]<<beginRc;
			endFiles[xs(hashBeginRc)%arg.nbFiles].write(reinterpret_cast<const char*>(&unitigIndice), 4);
		}

		end=unitig.substr(unitig.size()-arg.kmerSize,arg.kmerSize);
		endRc=revComp(end);
		hashEnd=str2num(end);
		hashEndRc=str2num(endRc);

		if(hashEnd<=hashEndRc){
			endFiles[xs(hashEnd)%arg.nbFiles]<<end;
			endFiles[xs(hashEnd)%arg.nbFiles].write(reinterpret_cast<const char*>(&unitigIndice), 4);
		}else{
			beginFiles[xs(hashEndRc)%arg.nbFiles]<<endRc;
			beginFiles[xs(hashEndRc)%arg.nbFiles].write(reinterpret_cast<const char*>(&unitigIndice), 4);
		}
		++unitigIndice;
	}
	if(arg.unitigThreshold==0){
		uint prev(0);
		for(uint i(1);i< abundance_unitigs.size(); ++i){
			cout<<abundance_unitigs[i]<<endl;
			if(prev!=0){
				if(abundance_unitigs[i]>=prev){
					arg.unitigThreshold=i-1;
					cout<<"	Unitig threshold chosed: "<<arg.unitigThreshold<<" "<<endl;
					break;
				}
				prev=abundance_unitigs[i];
			}else{
				prev=abundance_unitigs[i];
			}
		}
		for(uint i(0);i<unitigs.size();++i){
			if(coverages[i]<arg.unitigThreshold){
				unitigs[i]={};
			}else{
				break;
			}
		}
	}


	vector<bool> isaTip(unitigs.size(),false);


	if(arg.unitigThreshold==1){
		cout<<"\tTipping"<<endl;
		//TIPPING
		//FOREACH FILE
		uint passNumber(2);
		for(uint iPass(0);iPass<passNumber;++iPass){
			#pragma omp parallel for num_threads(arg.coreUsed)
			for(uint i=0; i< arg.nbFiles; ++i){
				string content,wordSeq,wordInt,seqBegin,seqEnd;
				vector<pair<string,uint32_t>> beginVector,endVector;
				uint32_t numberRead;
				uint sizeFileBegin(beginFiles[i].tellg());
				beginFiles[i].seekg(0,ios::beg);
				content.resize(sizeFileBegin);
				beginFiles[i].read(&content[0], content.size());
				for(uint ii(0); ii<content.size(); ){
					wordSeq=content.substr(ii,arg.kmerSize);
					wordInt=content.substr(ii+arg.kmerSize,4);
					memcpy(&numberRead,wordInt.c_str(),4);
					beginVector.push_back({wordSeq,numberRead});
					ii+=arg.kmerSize+4;
				}

				uint sizeFileEnd(endFiles[i].tellg());
				endFiles[i].seekg(0,ios::beg);
				content.resize(sizeFileEnd);
				endFiles[i].read(&content[0], content.size());
				for(uint ii(0); ii<content.size(); ){
					wordSeq=content.substr(ii,arg.kmerSize);
					wordInt=content.substr(ii+arg.kmerSize,4);
					memcpy(&numberRead,wordInt.c_str(),4);
					endVector.push_back({wordSeq,numberRead});
					ii+=arg.kmerSize+4;
				}
				sort(beginVector.begin(), beginVector.end());
				sort(endVector.begin(), endVector.end());
				uint indiceBegin(0), indiceEnd(0);
				vector<pair<uint,uint>> coverageComparison;
				while( not(indiceBegin==beginVector.size() and indiceEnd==endVector.size())){
					if(indiceBegin==beginVector.size()){
						seqBegin="Z";
					}else{
						seqBegin=beginVector[indiceBegin].first;
					}
					if(indiceEnd==endVector.size()){
						seqEnd="Z";
					}else{
						seqEnd=endVector[indiceEnd].first;
					}
					if(seqBegin==seqEnd){

						bool potential_crush(false);
						coverageComparison={};
						//if two begin and one is ten time larger remove the lower
						while(seqBegin==beginVector[indiceBegin].first){
							if(arg.ratioCoverage>0){
								if(unitigs[beginVector[indiceBegin].second].size()!=0){
									coverageComparison.push_back({coverages[beginVector[indiceBegin].second], abs(beginVector[indiceBegin].second)});
									if(coverages[beginVector[indiceBegin].second]<arg.high_coverage){
										potential_crush=true;
									}
								}
							}
							++indiceBegin;
							if(indiceBegin==beginVector.size()){
								break;
							}
						}

						if(potential_crush){

							for(uint iComp(0);iComp<coverageComparison.size();++iComp){
								for(uint iComp2(0);iComp2<coverageComparison.size();++iComp2){
									if(iComp!=iComp2  and not unitigs[coverageComparison[iComp].second].empty() and not unitigs[coverageComparison[iComp2].second].empty()){
										if(arg.ratioCoverage*coverageComparison[iComp].first < coverageComparison[iComp2].first and coverageComparison[iComp].first < arg.high_coverage){
											if(low_sistance_tip(unitigs[coverageComparison[iComp].second],unitigs[coverageComparison[iComp2].second],seqBegin)){
												unitigs[coverageComparison[iComp].second]={};
												advanced_tipping++;
											}
										}
									}
								}
							}




							//~ sort(coverageComparison.begin(),coverageComparison.end());
							//~ auto max_coverage(coverageComparison[coverageComparison.size()-1]);
							//~ for(uint iComp(0);iComp<coverageComparison.size()-1;++iComp){
								//~ if(passNumber==2){
									//~ if(isaTip[coverageComparison[iComp].second]){
										//~ if(low_sistance_tip(unitigs[coverageComparison[iComp].second],unitigs[max_coverage.second],seqBegin)){
											//~ #pragma omp critical(dataupdate)
											//~ {
												//~ unitigs[coverageComparison[iComp].second]={};
												//~ advanced_tipping++;
											//~ }
										//~ }
										//~ continue;
									//~ }
								//~ }
								//~ // LOW RELATIVE COVERAGE
								//~ if(arg.ratioCoverage*coverageComparison[iComp].first<max_coverage.first  and coverageComparison[iComp].first<arg.high_coverage){
									//~ //HAS THE SAME BEGIN/END THAN the BIG ONE
									//~ if(parallel_unitigs(unitigs[max_coverage.second],unitigs[coverageComparison[iComp].second],arg.kmerSize)){
										//~ #pragma omp critical(dataupdate)
										//~ {
											//~ unitigs[coverageComparison[iComp].second]={};
											//~ bulles++;
										//~ }
									//~ }else{
									//~ }
								//~ }
							//~ }
						}

						potential_crush=false;
						coverageComparison={};
						while(seqEnd==endVector[indiceEnd].first){
							if(arg.ratioCoverage>0){
								if(unitigs[endVector[indiceEnd].second].size()!=0){
									coverageComparison.push_back({coverages[endVector[indiceEnd].second],abs(endVector[indiceEnd].second)});
									if(coverages[endVector[indiceEnd].second]<arg.high_coverage){
										potential_crush=true;
									}
								}
							}
							++indiceEnd;
							if(indiceEnd==endVector.size()){
								break;
							}
						}

						if(potential_crush){

							for(uint iComp(0);iComp<coverageComparison.size();++iComp){
								for(uint iComp2(0);iComp2<coverageComparison.size();++iComp2){
									if(iComp!=iComp2  and not unitigs[coverageComparison[iComp].second].empty() and not unitigs[coverageComparison[iComp2].second].empty())
									if(arg.ratioCoverage*coverageComparison[iComp].first<coverageComparison[iComp2].first and coverageComparison[iComp].first<arg.high_coverage){
										if(low_sistance_tip(unitigs[coverageComparison[iComp].second],unitigs[coverageComparison[iComp2].second],seqEnd)){
											unitigs[coverageComparison[iComp].second]={};
											advanced_tipping++;
										}
									}
								}
							}


							//~ sort(coverageComparison.begin(),coverageComparison.end());
							//~ auto max_coverage(coverageComparison[coverageComparison.size()-1]);
							//~ for(uint iComp(0);iComp<coverageComparison.size()-1;++iComp){
								//~ if(passNumber==2){
									//~ if(isaTip[coverageComparison[iComp].second]){
										//~ if(low_sistance_tip(unitigs[coverageComparison[iComp].second],unitigs[max_coverage.second],seqBegin)){
											//~ #pragma omp critical(dataupdate)
											//~ {
												//~ unitigs[coverageComparison[iComp].second]={};
												//~ advanced_tipping++;
											//~ }
										//~ }
										//~ continue;
									//~ }
								//~ }
								//~ // LOW RELATIVE COVERAGE
								//~ if(arg.ratioCoverage*coverageComparison[iComp].first<max_coverage.first){
									//~ //HAS THE SAME BEGIN/END THAN the BIG ONE
										//~ if(parallel_unitigs(unitigs[max_coverage.second],unitigs[coverageComparison[iComp].second],arg.kmerSize)){
										//~ #pragma omp critical(dataupdate)
										//~ {
											//~ unitigs[coverageComparison[iComp].second]={};
											//~ bulles++;
										//~ }
									//~ }else{
									//~ }
								//~ }
							//~ }
						}
					}
					if(seqBegin<seqEnd){
						while(seqBegin==beginVector[indiceBegin].first){
							//TIP
							if(unitigs[beginVector[indiceBegin].second].size()<arg.tipingSize and unitigs[beginVector[indiceBegin].second].size()>0){
								//TIP WITH LOW COVERAGE
								if(coverages[beginVector[indiceBegin].second]<arg.low_coverage){
									#pragma omp critical(dataupdate)
									{
										++tiping;
										unitigs[beginVector[indiceBegin].second]={};
									}
								}else{
									//TIP BELOW THE SAVING THESHOLD
									if(coverages[beginVector[indiceBegin].second]<arg.high_coverage){
										//TIP CHECK low distance
										#pragma omp critical(TL)
										{
											isaTip[beginVector[indiceBegin].second]=true;
										}
									}
								}
							}
							++indiceBegin;
							if(indiceBegin==beginVector.size()){
								break;
							}
						}
						continue;
					}
					if(seqBegin>seqEnd){
						while(seqEnd==endVector[indiceEnd].first){
							//TIP
							if(unitigs[endVector[indiceEnd].second].size()<arg.tipingSize and unitigs[endVector[indiceEnd].second].size()>0){
								//TIP WITH LOW COVERAGE
								if(coverages[endVector[indiceEnd].second]<arg.low_coverage){
									#pragma omp critical(dataupdate)
									{
										++tiping;
										unitigs[endVector[indiceEnd].second]={};
									}
								}else{
									//TIP BELOW THE SAVING THESHOLD
									if(coverages[endVector[indiceEnd].second]<arg.high_coverage){
										//TIP CHECK low distance
										#pragma omp critical(TL)
										{
											isaTip[endVector[indiceEnd].second]=true;
										}
									}
								}

							}
							++indiceEnd;
							if(indiceEnd==endVector.size()){
								break;
							}
						}
						continue;
					}
				}
			}
		}
	}



	cout<<"\tRecompaction"<<endl;
	//RECOMPACTION
	#pragma omp parallel for num_threads(1)
	for(uint i=0; i< arg.nbFiles; ++i){
		string content,wordSeq,wordInt,seqBegin,seqEnd,str1,str2,compacted;
		vector<pair<string,uint32_t>> beginVector,endVector;
		uint numberRead;
		uint sizeFileBegin(beginFiles[i].tellg());
		beginFiles[i].seekg(0,ios::beg);
		content.resize(sizeFileBegin);
		beginFiles[i].read(&content[0], content.size());
		for(uint ii(0); ii<content.size(); ){
			wordSeq=content.substr(ii,arg.kmerSize);
			wordInt=content.substr(ii+arg.kmerSize,4);
			memcpy(&numberRead,wordInt.c_str(),4);
			if(not unitigs[numberRead].empty()){
				beginVector.push_back({wordSeq,numberRead});
			}else{
			}
			ii+=arg.kmerSize+4;
		}
		uint sizeFileEnd(endFiles[i].tellg());
		endFiles[i].seekg(0,ios::beg);
		content.resize(sizeFileEnd);
		endFiles[i].read(&content[0], content.size());
		for(uint ii(0); ii<content.size(); ){
			wordSeq=content.substr(ii,arg.kmerSize);
			wordInt=content.substr(ii+arg.kmerSize,4);
			memcpy(&numberRead,wordInt.c_str(),4);
			if(not unitigs[numberRead].empty()){
				endVector.push_back({wordSeq,numberRead});
			}
			ii+=arg.kmerSize+4;
		}
		sort(beginVector.begin(), beginVector.end());
		beginVector=uniqueOnly(beginVector);

		sort(endVector.begin(), endVector.end());
		endVector=uniqueOnly(endVector);
		uint indiceBegin(0), indiceEnd(0);
		while(indiceBegin<beginVector.size() and indiceEnd<endVector.size()){
			seqBegin=beginVector[indiceBegin].first;
			seqEnd=endVector[indiceEnd].first;
			if(seqBegin==seqEnd){
				#pragma omp critical(dataupdate)
				{
					uint position1(getPosition(unitigs,beginVector[indiceBegin].second));
					uint position2(getPosition(unitigs,endVector[indiceEnd].second));
					if(position1!=position2){
						str1=(bool2str(unitigs[position1]));
						str2=(bool2str(unitigs[position2]));
						compacted=(compaction(str1,str2,arg.kmerSize));
						unitigs[position1]=str2bool(compacted);
						coverages[position1]=((coverages[position1]*(str1.size()-arg.kmerSize))+(coverages[position2]*(str2.size()-arg.kmerSize)))/(str1.size()+str2.size()-2*arg.kmerSize);
						unitigs[position2]=int2bool(position1);
						compactions++;
					}
				}
				++indiceBegin;
				++indiceEnd;
				continue;
			}
			if(seqBegin<seqEnd){
				++indiceBegin;
				continue;
			}
			if(seqBegin>seqEnd){
				++indiceEnd;
				continue;
			}
		}
		remove((".begin"+to_string(i)).c_str());
		remove((".end"+to_string(i)).c_str());
	}

	//OUTPUT
	for(uint i(0); i<unitigs.size(); ++i){
		if((not unitigs[i].empty()) and (unitigs[i].size()%2==0) and coverages[i]>=arg.unitigThreshold){
			out<<">km:f:"<<coverages[i]<<"\n";
			out<<bool2str(unitigs[i])<<'\n';
		}else{
			unitigFiltered++;
		}
	}
	cout<<"\t..."<<endl;
	if(arg.unitigThreshold>1){
		cout<<"\tUnitig filtered: "+intToString(unitigFiltered)<<endl;
	}
	cout<<"\tTips removed: "+intToString(tiping)<<endl;
	cout<<"\tADVANCED Tips removed: "+intToString(advanced_tipping)<<endl;
	cout<<"\tIslands removed: "+intToString(island)<<endl;
	cout<<"\tBulle removed: "+intToString(bulles)<<endl;
	cout<<"\tUnitigs compacted: "+intToString(compactions)<<endl;

	auto endTime=system_clock::now();
    auto waitedFor=endTime-start;
    cout<<"\nCleaned in "<<duration_cast<seconds>(waitedFor).count()<<" seconds"<<endl;
    return (unitigFiltered+tiping+advanced_tipping+island+bulles);
}



//TODO multiple functions
//ONE FUNCTION TO RULE THEM ALL
int main(int argc, char ** argv){
	//INIT
	if(argc<3){
		usage();
		exit(0);
	}
	args_btt argument;
	argument.unitigThreshold=(1);
	string inputUnitig;
	//~ ifstream inUnitigs(inputUnitig);
	argument.tipingSize=(1);
	argument.kmerSize=(0);
	uint tipingStep(3);
	argument.coreUsed=(1);
	uint hashSize(8);//256 FILES
	argument.low_coverage=5;
	argument.ratioCoverage=(10);
	argument.unitigThreshold_MAX=(0);
	argument.high_coverage=10;
	string outFile("out_tipped");
	char c;
	while ((c = getopt (argc, argv, "u:k:t:c:h:f:a:o:T:m:l:L:")) != -1){
		switch(c){
		case 'u':
			inputUnitig=optarg;
			break;
		case 'k':
			argument.kmerSize=stoi(optarg);
			--argument.kmerSize;
			break;
		case 't':
			argument.tipingSize=2*stoi(optarg);
			break;
		case 'T':
			tipingStep=stoi(optarg);
			break;
		case 'c':
			argument.coreUsed=stoi(optarg);
			break;
		case 'h':
			hashSize=stoi(optarg);
			break;
		case 'f':
			argument.unitigThreshold=stoi(optarg);
			break;
		case 'm':
			argument.unitigThreshold_MAX=stoi(optarg);
			break;
		case 'o':
			outFile=optarg;
			break;
		case 'a':
			argument.ratioCoverage=stoi(optarg);
			break;
		case 'l':
			argument.low_coverage=stoi(optarg);
			break;
		case 'L':
			argument.high_coverage=stoi(optarg);
			break;
		}
	}
	argument.nbFiles=(1<<(hashSize-1));
	if(argument.kmerSize==0){
		cout<<"Please specify a kmer size"<<endl;
		exit(0);
	}

	for(uint i(1);i<=tipingStep;++i){
		cout<<"Step "<<i<<endl;
		if(i==tipingStep){
			int res(cleaning( outFile,  inputUnitig, argument));
			if(i>1){
				remove((outFile+to_string(i-1)).c_str());

			}
		}else{
			int res(cleaning( outFile+to_string(i),  inputUnitig, argument));
			if(res==0 and i>1){
				cout<<"nothing left to do I quit"<<endl;
				rename((outFile+to_string(i)).c_str(),outFile.c_str());
				return 0;
			}
			inputUnitig=outFile+to_string(i);
			if(i>1){
				remove((outFile+to_string(i-1)).c_str());
			}
		}
		argument.unitigThreshold=1;
	}

}
