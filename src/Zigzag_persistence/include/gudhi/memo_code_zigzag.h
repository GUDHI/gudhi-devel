
Implement the homology basis

Implement Morse complex

Implement general coefficients

//todo turn the public access, that is everywhere, into private.

//todo Documentation of compute_zigzag_persistence + syntax like compute_persistence in GUDHI.

//use a pool (boost?) ?


//FOR A POOL:
//In Zigzag_persistence creator:
//  , cell_pool_(new boost::object_pool< Cell > ()) //must be deleted
//As member of the class:
// boost::object_pool< Cell >                           * cell_pool_;
//In the code:
//        Cell * new_cell = new Cell(it2->key(), self1);    becomes,
//        Cell * new_cell = cell_pool_->construct(Cell(it2->key(), self1));
// and, for
//        Cell * tmp_ptr = &(*tmp_it); 
//        ...
//        delete tmp_ptr;                                  becomes,
//        cell_pool_->free(tmp_ptr);




//todo interval compare and interval display policies (see Persistent_cohomology.h)

//todo make all fields private, watch for std::swap functions





      // //add all chains with smaller <d death and larger <b birth than max_avail_b
      // for(auto chain_passed_it = chains_in_F.rbegin();//all with smaller <d death
      //          chain_passed_it != chain_f_it;  ++chain_passed_it) 
      // {//but
      //   if(cmp_birth((*chain_passed_it)->birth(), max_avail_b)) 
      //   {//larger <b birth
      //     plus_equal_column( (*chain_f_it), (*chain_f_it)->column(), (*chain_passed_it)->column() );
      //   } 
      // }

Persistence diagram:

end of function backward arrow, case remove someone in H, G becomes F
    else { //in H    -> paired with c_g, that now belongs to F now
      curr_col->paired_col_->assign_paired_chain(nullptr);
      curr_col->paired_col_->assign_birth(num_arrow_); //closed interval
    }





// Simplex_tree.h

//   Cofaces_simplex_range cofaces_simplex_range(const Simplex_handle simplex, int codimension) {
//     return cofaces_simplex_range(simplex, codimension,
//                               // std::integral_constant<bool, Options::link_nodes_by_label>{}
//                               std::integral_constant<bool, false>{}
//                               );
//   }



  Cofaces_simplex_range cofaces_simplex_range(const Simplex_handle simplex, int codimension, std::false_type) {
    Cofaces_simplex_range cofaces;
    // codimension must be positive or null integer
    assert(codimension >= 0);
    Simplex_vertex_range rg = simplex_vertex_range(simplex);
    std::vector<Vertex_handle> copy(rg.begin(), rg.end());
    if (codimension + static_cast<int>(copy.size()) > dimension_ + 1 ||
        (codimension == 0 && static_cast<int>(copy.size()) > dimension_))  // n+codimension greater than dimension_
      return cofaces;
    // must be sorted in decreasing order
    assert(std::is_sorted(copy.begin(), copy.end(), std::greater<Vertex_handle>()));
    bool star = codimension == 0;

    rec_coface(copy, &root_, 1, cofaces, star, codimension + static_cast<int>(copy.size()));
    return cofaces;
  }

































//former arrow transposition case stuy, with chains in G


//     switch( curr_col->birth() ) 
//     {
//       case -2://====================================================================
//       {//case H x *       c_s is in H
//         switch( other_col->birth() ) 
//         { 
//           case -2://----------------------------------------------------------------
//           {//Case H x H, c_s+c_t paired w/ c_g+c_g', of death max<d{g,g'}=max
//             auto curr_p_col  = curr_col->paired_col_; //c_s paired with c_g
//             auto other_p_col = other_col->paired_col_;//c_t paired with c_g'
//             if( curr_p_col->lowest_idx_ < other_p_col->lowest_idx_) 
//             {//g<g', -->|c_s|
//               plus_equal_column( other_p_col
//                                , other_p_col->column()//c_g' <- c_g+c_g', of death
//                                , curr_p_col->column() );// g', low idx same, etc
//               plus_equal_column( other_col
//                                 , other_col->column()
//                                , curr_col->column() ); //c_t <- c_t+c_s
                
//                // ++_other_pcol_is_c1_;
//                // ++_other_col_is_c1_;

//                          #ifdef _PROFILING_ZIGZAG_PERSISTENCE_ 
//                            ++_num_arrow_trans_case_study_plus_equal_col_;
//                            ++_num_arrow_trans_case_study_plus_equal_col_;
//                            ++case1;
//                          #endif

//                return curr_col;//continue with c_s, paired with c_g of min death g
//                }
//                else 
//                {// g' < g, continue with --> |c_t|
//                plus_equal_column( curr_p_col
//                                 , curr_p_col->column()    //c_g <- c_g+c_g', death
//                                 , other_p_col->column() );// g, low idx same, etc
//                plus_equal_column( curr_col, curr_col->column()//c_s <- c_s+c_t, of
//                                 , other_col->column());       //low idx t now
//                //exchange lowest_idx, update lowidx_to_matidx structure
//                exchange_lowest_indices_chains(curr_col, other_col);
                  
//                // ++_curr_pcol_is_c1_;
//                // ++_curr_col_is_c1_;
//                          #ifdef _PROFILING_ZIGZAG_PERSISTENCE_ 
//                            ++_num_arrow_trans_case_study_plus_equal_col_;
//                            ++_num_arrow_trans_case_study_plus_equal_col_;
//                            ++case2;
//                          #endif

//                return other_col; //continue with c_t, paired w. c_g' of min death g'
//                }
//          }//end case H x H
//           default: 
//           { //in H x (F U G) : c_s+c_t in H, paired with c_g
//              // auto curr_p_col = curr_col->paired_col_; //c_s paired with c_g
//                plus_equal_column(curr_col, curr_col->column() //(still in H) <-> c_g
//                                 , other_col->column());       //c_s <- c_s+c_t
//                //exchange lowest_idx, update lowidx_to_matidx structure
//                exchange_lowest_indices_chains(curr_col, other_col);
//                //exchange all members, except the column_ pointers, connected to cells
//                // exchange_pairings_chains(curr_col, other_col);
                  
//                // ++_curr_col_is_c1_;
//                      #ifdef _PROFILING_ZIGZAG_PERSISTENCE_ 
//                        ++_num_arrow_trans_case_study_plus_equal_col_;
//                        ++case3;
//                      #endif

//                return other_col; //continue with c_t, still in F
//          }//end case H x (F U G)
//         }//end nested switch
//       }//end case H x *
   
      

//       case -1://====================================================================
//       {//case G x *       c_s is in G
//         switch( other_col->birth() ) 
//         {
//           case -2://----------------------------------------------------------------
//           { //in G x H
//                plus_equal_column( other_col, other_col->column()
//                                 , curr_col->column() );

//                      #ifdef _PROFILING_ZIGZAG_PERSISTENCE_ 
//                        ++_num_arrow_trans_case_study_plus_equal_col_;
//                        ++case4;
//                      #endif

//                return curr_col;
//         }//end case G x H
//         case -1://------------------------------------------------------------------
//         {//in G x G,          c_s+c_t (in G) has death max<d{h,h'}=max{h,h'}
//            auto curr_p_col  = curr_col->paired_col_;  //c_s paired with c_h
//            auto other_p_col = other_col->paired_col_; //c_t paired with c_h'
//            if( curr_p_col->lowest_idx_ < other_p_col->lowest_idx_ ) 
//            {// h < h'
//            plus_equal_column( other_p_col, other_p_col->column() //c_h' <- c_h+c_h'
//                             , curr_p_col->column() );     //of death h', low idx h'
//            plus_equal_column( other_col, other_col->column()//c_t <- c_s+c_t
//                             , curr_col->column() );         //paired with c_h'

//                      #ifdef _PROFILING_ZIGZAG_PERSISTENCE_ 
//                        ++_num_arrow_trans_case_study_plus_equal_col_;
//                        ++_num_arrow_trans_case_study_plus_equal_col_;
//                        ++case5;
//                      #endif
//             return curr_col;//continue with c_s, still paired with c_h
//            }//endif h < h'
//            else 
//            {//h' < h
//            plus_equal_column( curr_p_col, curr_p_col->column()//c_h <- c_h+c_h'
//                             , other_p_col->column() );//of death h, low idx h
//            plus_equal_column( curr_col, curr_col->column()//still paired with c_h
//                             , other_col->column());//c_s <- c_s+c_t, of low idx t
//            //exchange lowest_idx, update lowidx_to_matidx structure
//            exchange_lowest_indices_chains(curr_col, other_col);

//                      #ifdef _PROFILING_ZIGZAG_PERSISTENCE_ 
//                        ++_num_arrow_trans_case_study_plus_equal_col_;
//                        ++_num_arrow_trans_case_study_plus_equal_col_;
//                        ++case6;
//                      #endif
//            return other_col; //continue with c_t, of min death h' and low idx s
//            }//end else h' < h
//         }//end case G x G
//         default://-----------------------------------------------------------------
//         { //in G x F, c_s+c_t (in F) has maximal birth
//            plus_equal_column( other_col, other_col->column()
//                             , curr_col->column() );//c_t <- c_t+c_s, in F

//                  #ifdef _PROFILING_ZIGZAG_PERSISTENCE_       
//                    ++_num_arrow_trans_case_study_plus_equal_col_;
//                    ++case7;
//                  #endif

//            return curr_col;
//         }//end case GF
//       }//end nested switch
//     }//end case G x *

//     default://======================================================================
//     {//case F x *         c_s is in F
//       switch( other_col->birth() ) 
//       {
//        case -2://------------------------------------------------------------------
//        { // in F x H
//            plus_equal_column( other_col, other_col->column()
//                               , curr_col->column() );//c_t <- c_s+c_t

//                  #ifdef _PROFILING_ZIGZAG_PERSISTENCE_ 
//                    ++_num_arrow_trans_case_study_plus_equal_col_;
//                    ++case8;
//                  #endif

//            return curr_col;
//         }//end case F x H
//         case -1://------------------------------------------------------------------
//         { //in F x G:       c_s+c_t is in F, of max birth
//            plus_equal_column(curr_col, curr_col->column()//c_s <- c_s+c_t, still
//                              , other_col->column());// in F, low idx t
//            //exchange lowest_idx, update lowidx_to_matidx structure
//            exchange_lowest_indices_chains(curr_col, other_col);
      
//                  #ifdef _PROFILING_ZIGZAG_PERSISTENCE_ 
//                    ++_num_arrow_trans_case_study_plus_equal_col_;
//                    ++case9;
//                  #endif

//           return other_col; //still in G
//         }// end case F x G
//         default://------------------------------------------------------------------
//         { //in F x F:       c_s+c_t has max<=b birth
//          if(birth_ordering_.birth_order(curr_col->birth(), other_col->birth())) 
//          {
//             plus_equal_column( other_col, other_col->column()     //b_s <b b_t
//                                     , curr_col->column() );//c_t <- c_s+c_t, of low idx t
        
//                      #ifdef _PROFILING_ZIGZAG_PERSISTENCE_ 
//                        ++_num_arrow_trans_case_study_plus_equal_col_;
//                        ++case10;
//                      #endif

//            return curr_col;
//           }//endif
//           else 
//           { //b_t <b b_s
//              plus_equal_column(curr_col, curr_col->column()
//                               , other_col->column());//c_s <- c_s+c_t, of low idx t
//              //exchange lowest_idx, update lowidx_to_matidx structure
//              exchange_lowest_indices_chains(curr_col, other_col);

//                      #ifdef _PROFILING_ZIGZAG_PERSISTENCE_ 
//                        ++_num_arrow_trans_case_study_plus_equal_col_;
//                        ++case11;
//                      #endif

//              return other_col;
//            }//end else
//         }//end case F x F
//       }//end nested switch
//     }//end case F x *
//   }//end switch
// }//end function
















    //persistence diagram
      //   int curr_dim = -2;//tmp_diag.begin()->dim_;
  //   // int curr_num_intervals = num_intervals;

  //   int i;//=0;
  //   os << "#i  birth  death     length \n";
  //   for(auto bar : tmp_diag) {
  //     if(curr_dim != bar.dim_) {
  //       // std::cout << "------------------------------- dim " << bar.dim_ << " \n";
  //       // curr_num_intervals = num_intervals; 
  //       curr_dim = bar.dim_; 
  //       i=0;
  //     }
  //    //  if(curr_num_intervals > 0 && bar.length() > 0.0001) {
  //     // --curr_num_intervals;
  //     if(bar.b_ < bar.d_) {
  //       // os << bar.dim_ << "   " << bar.b_ << " " << bar.d_ <<
  //       //  "      " << bar.length() << "        #" << i++ << " \n";
  //       os << "#" << i++ << "  " << bar.b_ << " " << bar.d_ <<
  //        "             " << bar.length() << " \n";
  //     }
  //     else {
  //      //  os << bar.dim_ << "   " << bar.d_ << " " << bar.b_ <<
  //        // "      " << bar.length()  << "        #" << i++ << " \n";
  //       os << "#" << i++ << "  " << bar.d_ << " " << bar.b_ <<
  //        "             " << bar.length() << " \n";
  //     }
  //   }
  // }
  





  // //Dionysus: In interval
  // // (idx i)_R_{eta eps_i}(Pi) -> .. \leftarrow R_{eta eps_i+1}(Pi+1)_(idx j), anything
  // // created from index i+1 on has birth eps_i+1, anything destroyed has death eps_i
  // void output_diagram_log2( std::ostream& os = std::cout
  //                         , Filtration_value shortest_interval = 0.) {

  //   std::stable_sort(filtration_values_.begin(), filtration_values_.end(),
  //        []( std::pair< Simplex_key, Filtration_value > p1
  //          , std::pair< Simplex_key, Filtration_value > p2 )
  //          { return p1.first < p2.first; }
  //   );
  //   std::vector< interval_t > tmp_diag;
  //   tmp_diag.reserve(persistence_diagram_.size());
  //   for(auto bar : persistence_diagram_)
  //   {
  //     auto it_b = //lower_bound(x) returns leftmost y s.t. x <= y
  //       std::lower_bound( filtration_values_.begin(), filtration_values_.end()
  //             , std::pair<Simplex_key, Filtration_value>(bar.b_
  //                                       , std::numeric_limits<double>::infinity() )
  //             , []( std::pair<Simplex_key, Filtration_value> p1
  //                 , std::pair<Simplex_key, Filtration_value> p2) 
  //                 {
  //                   return p1.first < p2.first; 
  //                 }
  //       );
  //     //
  //     if(it_b == filtration_values_.end() || it_b->first > bar.b_) { --it_b; }

  //     auto it_d = //upper_bound(x) returns leftmost y s.t. x < y, or last
  //       std::upper_bound( filtration_values_.begin(), filtration_values_.end()
  //             , std::pair<Simplex_key, Filtration_value>(bar.d_
  //                                        , std::numeric_limits<double>::infinity() )
  //             , []( std::pair<Simplex_key, Filtration_value> p1
  //                 , std::pair<Simplex_key, Filtration_value> p2) 
  //                 {
  //                   return p1.first < p2.first; 
  //                 }
  //       );
  //     //discard interval strictly included between two consecutive indices
  //     --it_d;
  //     // if(it_b->second != it_d->second) {
  //     if( std::abs(log2(it_b->second) - log2(it_d->second)) >= shortest_interval ) {
  //       tmp_diag.emplace_back(bar.dim_, log2(it_b->second), log2(it_d->second) );
  //     }
  //   }
  //   cmp_intervals_by_length cmp;
  //   std::stable_sort(tmp_diag.begin(), tmp_diag.end(), cmp);

  //   if(tmp_diag.empty()) {return;}

  //   int curr_dim = -2;//tmp_diag.begin()->dim_;
  //   // int curr_num_intervals = num_intervals;

  //   int i;//=0;
  //   os << "#i  birth  death     length \n";
  //   for(auto bar : tmp_diag) {
  //     if(curr_dim != bar.dim_) {
  //       // std::cout << "------------------------------- dim " << bar.dim_ << " \n";
  //       // curr_num_intervals = num_intervals; 
  //       curr_dim = bar.dim_; 
  //       i=0;
  //     }
  //    //  if(curr_num_intervals > 0 && bar.length() > 0.0001) {
  //     // --curr_num_intervals;
  //     if(bar.b_ < bar.d_) {
  //       // os << bar.dim_ << "   " << bar.b_ << " " << bar.d_ <<
  //       //  "      " << bar.length() << "        #" << i++ << " \n";
  //       os << "#" << i++ << "  " << bar.b_ << " " << bar.d_ <<
  //        "             " << bar.length() << " \n";
  //     }
  //     else {
  //      //  os << bar.dim_ << "   " << bar.d_ << " " << bar.b_ <<
  //        // "      " << bar.length()  << "        #" << i++ << " \n";
  //       os << "#" << i++ << "  " << bar.d_ << " " << bar.b_ <<
  //        "             " << bar.length() << " \n";
  //     }
  //   }
  // }




  //   void stat_mat() {
  //     size_t size_mat = matrix_.size();
  //     std::cout << "      *** " << " #columns = " << size_mat << ",  ";

  //     std::vector<size_t> sizes;
  //     sizes.reserve(size_mat);

  //     std::vector<size_t> sizesF;
  //     sizesF.reserve(size_mat);
  //     std::vector<size_t> sizesG;
  //     sizesG.reserve(size_mat);      
  //     std::vector<size_t> sizesH;
  //     sizesH.reserve(size_mat);


  //     unsigned long long num_cells = 0;
  //     unsigned long long num_cellsF = 0;
  //     unsigned long long num_cellsG = 0;
  //     unsigned long long num_cellsH = 0;

  //     for(auto & col : matrix_) {
  //       size_t size_col = col.column().size();
  //       // size_t size_col = 0;
  //       // for(auto c : *(col.column())) { ++size_col; }
  //       num_cells += size_col;
  //       sizes.push_back(size_col);

  //       if(col.inF()) { sizesF.push_back(size_col); num_cellsF += size_col; }
  //       if(col.inG()) { sizesG.push_back(size_col); num_cellsG += size_col; }
  //       if(col.inH()) { sizesH.push_back(size_col); num_cellsH += size_col; }
  //     }

  //     std::cout << "#cellsTotal = " << num_cells << ",  ";
  //     std::cout << "av. col size = " << (double)num_cells/(double)size_mat << " ***\n";

  //     if(size_mat > 10) {
  //       std::stable_sort(sizes.begin(),sizes.end());
  //       std::cout << "  " << sizes.size() << "     *** " << "10%-" << sizes[(size_t)(size_mat/10)] << "  50%-" << sizes[(size_t)(size_mat/2)] << "  90%-" << sizes[(size_t)(9*size_mat/10)] << "  95%-" << sizes[(size_t)(95*size_mat/100)] << "  98%-" << sizes[(size_t)(98*size_mat/100)] << "  99%-" << sizes[(size_t)(99*size_mat/100)] << " \n";
  //     }

  //     std::cout << "#cellsF = " << num_cellsF << ",  ";
  //     std::cout << "av. colF size = " << (double)num_cellsF/(double)sizesF.size() << " ***\n";
  //     if(sizesF.size() > 10) {
  //       std::stable_sort(sizesF.begin(),sizesF.end());
  //       std::cout << "  " << sizesF.size() << "     *** " << "10%-" << sizesF[(size_t)(sizesF.size()/10)] << "  50%-" << sizesF[(size_t)(sizesF.size()/2)] << "  90%-" << sizesF[(size_t)(9*sizesF.size()/10)] << "  95%-" << sizesF[(size_t)(95*sizesF.size()/100)] << "  98%-" << sizesF[(size_t)(98*sizesF.size()/100)] << "  99%-" << sizesF[(size_t)(99*sizesF.size()/100)] << "  100%-" << sizesF[sizesF.size()-1] << " \n";
  //     }
  //     std::cout << "#cellsG = " << num_cellsG << ",  ";
  //     std::cout << "av. colG size = " << (double)num_cellsG/(double)sizesG.size() << " ***\n";
  //     if(sizesG.size() > 10) {
  //       std::stable_sort(sizesG.begin(),sizesG.end());
  //       std::cout << "  " << sizesG.size() << "     *** " << "10%-" << sizesG[(size_t)(sizesG.size()/10)] << "  50%-" << sizesG[(size_t)(sizesG.size()/2)] << "  90%-" << sizesG[(size_t)(9*sizesG.size()/10)] << "  95%-" << sizesG[(size_t)(95*sizesG.size()/100)] << "  98%-" << sizesG[(size_t)(98*sizesG.size()/100)] << "  99%-" << sizesG[(size_t)(99*sizesG.size()/100)] << "  100%-" << sizesG[sizesG.size()-1] << " \n";
  //     }
  //     std::cout << "#cellsH = " << num_cellsH << ",  ";
  //     std::cout << "av. colH size = " << (double)num_cellsH/(double)sizesH.size() << " ***\n";
  //     if(sizesH.size() > 10) {
  //       std::stable_sort(sizesH.begin(),sizesH.end());
  //       std::cout << "  " << sizesH.size() << "     *** " << "10%-" << sizesH[(size_t)(sizesH.size()/10)] << "  50%-" << sizesH[(size_t)(sizesH.size()/2)] << "  90%-" << sizesH[(size_t)(9*sizesH.size()/10)] << "  95%-" << sizesH[(size_t)(95*sizesH.size()/100)] << "  98%-" << sizesH[(size_t)(98*sizesH.size()/100)] << "  99%-" << sizesH[(size_t)(99*sizesH.size()/100)] << "  100%-" << sizesH[sizesH.size()-1] << " \n";
  //     }
  //   }

  //   void display_mat() {
  //     std::cout << "---------------------------beg \n";
  //     for(auto & col : matrix_)
  //     {
  //       // if(col.paired_col_ == nullptr && col.birth() != col.lowest_idx_) {
  //       //   std::cout << "birth != lowest_idx !!\n";
  //       // } <---- birth != lowest_idx
  //       std::cout << " : [low_idx=" << col.lowest_idx_ <<"] ";
  //       std::cout << "[birth=" << col.birth_ <<"] ";
  //       if(col.paired_col_ != nullptr) { 
  //         std::cout << "[paired with=" << col.paired_col_->lowest_idx_ <<"]     ";
  //       }
  //       else { std::cout << "[paired with=" << "F" <<"]     ";}
  //       std::cout <<"| ";
  //       for(auto &cell : *(col.column_)) { std::cout << cell.key() << " ; "; }
  //       std::cout << " |" << std::endl;
  //     }
  //     std::cout << "----------- \n";
  //     // std::cout << "Bv: "; for(auto b : birth_vector_.inv_b_vec) { std::cout << b << " "; } std::cout << "\n";
  //     // std::cout << "Bv: "; for(auto pp : birth_ordering_.birth_to_pos_) { std::cout << "("<<pp.first<<";"<<pp.second<<") ";}
  //     // std::cout << "    maxb = " << birth_ordering_.max_birth_pos_ << " , minb = " << birth_ordering_.min_birth_pos_ << "\n";
  //     std::cout << "Persistence diagram: "; 
  //     for(auto bd: persistence_diagram_)
  //     { 
  //       std::cout << "["<<bd.b_<<";"<<bd.d_<< "  d" << bd.dim_ << "] "; 
  //     } 
  //     std::cout << "\n";
  //     // std::cout << "l_to_m: "; for(auto ltm : lowidx_to_matidx_)
  //     //                       { std::cout << ltm.first << "->" << ltm.second->lowest_idx_ << " ";} std::cout << "\n";
  //     for(auto ltm : lowidx_to_matidx_) {
  //       if(ltm.first != ltm.second->lowest_idx_) {std::cout << "ERROR ltm != lowest_idx ??\n";}
  //     }
  //     std::cout << "---------------------------end \n";
  //   }








//OLD arrow_transposition

//return the new value of curr_col
// matrix_chain * arrow_transposition_case_study( matrix_chain * curr_col
//                                              , matrix_chain * other_col )
// {
//   switch( curr_col->birth() ) 
//   {
//     case -2: 
//     { //in H
//       switch( other_col->birth() ) 
//       {
//         case -2: {//Case H x H, c_s+c_t paired w. c_g+c_g', of death max<d{g,g'}=max
//           auto curr_p_col  = curr_col->paired_col_; //c_s paired with c_g
//           auto other_p_col = other_col->paired_col_;//c_t paired with c_g'
//           if( curr_p_col->lowest_idx_ < other_p_col->lowest_idx_) {//g<g', -->|c_s|
//             plus_equal_column( other_p_col, other_p_col->column() 
//                                           , curr_p_col->column() );//c_g' <- c_g+c_g'
//             plus_equal_column( other_col, other_col->column() 
//                                         , curr_col->column() ); //c_t <- c_t+c_s
//           }
//           else {// g' < g, continue with --> |c_t|
//             plus_equal_column( curr_p_col, curr_p_col->column() 
//                              , other_p_col->column() );//c_g <- c_g+c_g'
//             plus_equal_column( curr_col, curr_col->column()
//                                        , other_col->column());//c_s <- c_s+c_t      
//             return other_col; //continue with c_t


//             // plus_equal_column( curr_p_col, curr_p_col->column() 
//             //                              , other_p_col->column() );//c_g <- c_g+c_g'
//             // col_swap( curr_col , other_col ); //exhange columns so as h' survives
//             // //change pairing
//             // std::swap( curr_p_col->paired_col_, other_p_col->paired_col_ );
//             // std::swap( curr_col->paired_col_, other_col->paired_col_ );
//             // plus_equal_column( other_col, other_col->column() , curr_col->column() ); //h' <- h+h'
//           }
//           break;
//         }




//         case -1: { //in H x G, (in H) c_s+c_t <-> c_g 
//           auto curr_p_col = curr_col->paired_col_; //c_s paired with c_g
//           plus_equal_column(curr_col, curr_col->column() //(still in H) <-> c_g
//                                     , other_col->column());//c_s <- c_s+c_t 
//           std::swap(curr_col, other_col); //continue with c_t, still in G      
//           break;


//           // auto curr_p_col = curr_col->paired_col_; 
//           // auto other_p_col = other_col->paired_col_;

//           // col_swap( curr_col , other_col );
//           // std::swap( curr_p_col->paired_col_, other_p_col->paired_col_ );
//           // std::swap( curr_col->paired_col_, other_col->paired_col_ );

//           // plus_equal_column( other_col, other_col->column() , curr_col->column() );

//           // other_col->assign_birth(-2); curr_col->assign_birth(-1);

//           // break;
//         }
//         default: 
//         { //in H x F
//           auto curr_p_col = curr_col->paired_col_; //c_s paired with c_g
//           plus_equal_column(curr_col, curr_col->column() //(still in H) <-> c_g
//                                     , other_col->column());//c_s <- c_s+c_t 
//           std::swap(curr_col, other_col); //continue with c_t, still in F      
//           break;

//           // col_swap( curr_col , other_col );
//           // plus_equal_column( other_col, other_col->column() , curr_col->column() );
//           // std::swap( curr_col->paired_col_, other_col->paired_col_ );
//           // other_col->paired_col_->paired_col_ = other_col;
//           // std::swap( curr_col->birth_ , other_col->birth_ );
//           // break;
//         }
//       }
//       break;
//     }
//     case -1: 
//     { //in G


//  //     std::cout << "Transposition with column in G \n \n \n ";

//       switch( other_col->birth() ) 
//       {
//         case -2: { //in G x H ok

//           // std::cout << "GH ";
//           // if(other_col->paired_col_ == curr_col) 
//           //   { std::cout << " Complain.......................................\n";}

//           plus_equal_column( other_col, other_col->column() , curr_col->column() );
//           break;
//         }
//         case -1: 
//         { //in G x G 

//           auto curr_p_col = curr_col->paired_col_; auto other_p_col = other_col->paired_col_;
//           if( curr_p_col->lowest_idx_ < other_p_col->lowest_idx_ ) 
//           {
//             plus_equal_column( other_p_col, other_p_col->column() , curr_p_col->column() ); //h' <- h+h'
//             plus_equal_column( other_col, other_col->column() , curr_col->column() ); //g' <- g+g'
//           }
//           else 
//           {
//             plus_equal_column( curr_p_col, curr_p_col->column() , other_p_col->column() ); 
//             col_swap( curr_col , other_col ); 
//             //change pairing
//             std::swap( curr_p_col->paired_col_, other_p_col->paired_col_ );
//             std::swap( curr_col->paired_col_, other_col->paired_col_ );

//             plus_equal_column( other_col, other_col->column() , curr_col->column() ); 
//           }
//           break;
//         }
//         default: 
//         { //in G x F
//           //---leave c1+c2 (belongs to F), continue with c1 in G
//           plus_equal_column( other_col, other_col->column() , curr_col->column() );
//           break;
//         }
//         break;
//       }
//     }   
//     default: 
//     { // in F
//       switch( other_col->birth() ) 
//       {
//         case -2: { // in F x H


//           plus_equal_column( other_col, other_col->column() , curr_col->column() );                  
//           break;
//         }
//         case -1: 
//         { //in F x G:
//           //---leave c1+c2 (belongs to F), continue with c2 in G
//           col_swap( curr_col , other_col ); 
//           //change pairing
//           std::swap( curr_col->paired_col_, other_col->paired_col_ );
//           curr_col->paired_col_->paired_col_ = curr_col;
//           std::swap( curr_col->birth_ , other_col->birth_ );
//           break;
//         }
//         default: 
//         { //in F x F: 
//           //---leave c1+c2 (has max birth), continue with c1 or c2 with min birth

//           if(birth_ordering_.birth_order(curr_col->birth(), other_col->birth()))
//        // if( birth_vector_[curr_col->birth()] < birth_vector_[other_col->birth()] )
//           {plus_equal_column( other_col, other_col->column() , curr_col->column() );}
//           else 
//           { 
//             col_swap( curr_col , other_col ); 
//             std::swap( curr_col->birth_ , other_col->birth_ );
//             plus_equal_column( other_col, other_col->column() , curr_col->column() );
//           }
//           break;
//         }
//         break;
//       }
//     }  
//   } //end switch
// }




































// //return the new value of curr_col
// matrix_chain * arrow_transposition_case_study( matrix_chain * curr_col, matrix_chain * other_col )
// {
//   switch( curr_col->birth_ ) 
//   {
//     case -2: 
//     { //in H
//       switch( other_col->birth_ ) 
//       {
//         case -2: 
//         { //in H x H
//           if( curr_col->paired_col_->lowest_idx_ < other_col->paired_col_->lowest_idx_) 
//           { //g < g'
//             plus_equal_column( other_col->paired_col_->column()
//                              , curr_col->paired_col_->column()); //g' <- g+g'
//             plus_equal_column( other_col->column()
//                              , curr_col_->column() ); //h' <- h+h'


//             matrix_[matrix_[curr_col].paired_col].paired_col = curr_col+1;//new pairing gs
//             matrix_[matrix_[curr_col+1].paired_col].paired_col = curr_col;
//             std::swap(matrix_[curr_col], matrix_[curr_col+1]);          //permute columns
//             plus_equal_column(matrix_[curr_col].column, matrix_[curr_col+1].column); //h+h'
//           }
//           else {
//             plus_equal_column( matrix_[matrix_[curr_col].paired_col].column
//                              , matrix_[matrix_[curr_col+1].paired_col].column ); //g+g'
//             plus_equal_column( matrix_[curr_col].column
//                              , matrix_[curr_col+1].column );//h+h'
//             std::swap(matrix_[curr_col].lowest_idx, matrix_[curr_col+1].lowest_idx);
//             std::swap(matrix_[curr_col].row, matrix_[curr_col+1].row);
//           }
//           break;
//         }
//         case -1: { //in H x G
//           plus_equal_column( matrix_[curr_col].column
//                            , matrix_[curr_col+1].column );//g+h'
//           std::swap(matrix_[curr_col].lowest_idx, matrix_[curr_col+1].lowest_idx);
//           std::swap(matrix_[curr_col].row, matrix_[curr_col+1].row);
//           break;
//         }
//         default: { //in H x F
//           plus_equal_column( matrix_[curr_col].column
//                            , matrix_[curr_col+1].column );//f+h
//           std::swap(matrix_[curr_col].lowest_idx, matrix_[curr_col+1].lowest_idx);
//           std::swap(matrix_[curr_col].row, matrix_[curr_col+1].row);
//           break;
//         }
//       }
//       break;
//     }
//     case -1: 
//     { //in G
//       switch(matrix_[curr_col+1].birth) 
//       {
//         case -2: { //in G x H
//           matrix_[matrix_[curr_col].paired_col].paired_col = curr_col+1; //new pairing
//           matrix_[matrix_[curr_col+1].paired_col].paired_col = curr_col;
//           std::swap(matrix_[curr_col], matrix_[curr_col+1]);          //permute columns
//           plus_equal_column( matrix_[curr_col].column
//                            , matrix_[curr_col+1].column );//h'+g                  
//           break;
//         }
//         case -1: 
//         { //in G x G
//           if(matrix_[curr_col].paired_col < matrix_[curr_col+1].paired_col) {
//             plus_equal_column(  matrix_[matrix_[curr_col+1].paired_col].column
//                               , matrix_[matrix_[curr_col].paired_col].column ); //h+h'
//             matrix_[matrix_[curr_col].paired_col].paired_col = curr_col+1;//new pairing hs
//             matrix_[matrix_[curr_col+1].paired_col].paired_col = curr_col;
//             std::swap(matrix_[curr_col], matrix_[curr_col+1]);          //permute columns
//             plus_equal_column(matrix_[curr_col].column, matrix_[curr_col+1].column); //g+g'
//           }
//           else {
//             plus_equal_column( matrix_[matrix_[curr_col].paired_col].column
//                              , matrix_[matrix_[curr_col+1].paired_col].column ); //h+h'
//             plus_equal_column( matrix_[curr_col].column
//                              , matrix_[curr_col+1].column );//g+g'
//             std::swap(matrix_[curr_col].lowest_idx, matrix_[curr_col+1].lowest_idx);
//             std::swap(matrix_[curr_col].row, matrix_[curr_col+1].row);
//           }
//           break;
//         }
//         default: 
//         { //in G x F
//           matrix_[matrix_[curr_col].paired_col].paired_col = curr_col+1;//new pairing hs
//           std::swap(matrix_[curr_col], matrix_[curr_col+1]);          //permute columns
//           plus_equal_column( matrix_[curr_col].column
//                            , matrix_[curr_col+1].column );//f+g
//           break;
//         }
//         break;
//       }
//     }   
//     default: 
//     { // in F
//       switch(matrix_[curr_col+1].birth) 
//       {
//         case -2: { // in F x H
//           matrix_[matrix_[curr_col+1].paired_col].paired_col = curr_col; //new pairing
//           std::swap(matrix_[curr_col], matrix_[curr_col+1]);          //permute columns
//           plus_equal_column( matrix_[curr_col].column
//                            , matrix_[curr_col+1].column );//h+f                  
//           break;
//         }
//         case -1: { //in F x G
//           plus_equal_column( matrix_[curr_col].column
//                            , matrix_[curr_col+1].column );//f+g
//           std::swap(matrix_[curr_col].lowest_idx, matrix_[curr_col+1].lowest_idx);
//           std::swap(matrix_[curr_col].row, matrix_[curr_col+1].row);
//           break;
//         }
//         default: 
//         { //in F x F
//           if( birth_vector_[matrix_[curr_col].birth] < birth_vector_[matrix_[curr_col+1].birth] ) 
//           { std::swap(matrix_[curr_col], matrix_[curr_col+1]); }  //permute columns
//           else 
//             { std::swap(matrix_[curr_col].lowest_idx, matrix_[curr_col+1].lowest_idx); 
//               std::swap(matrix_[curr_col].row, matrix_[curr_col+1].row); 
//             }
//           plus_equal_column( matrix_[curr_col].column
//                            , matrix_[curr_col+1].column );//f+f'
//           break;
//         }
//         break;
//       }
//     }  
//   } //end switch
// }








//   std::map< Simplex_key, Column * >          mat_; // matrix

//   // to update properly
//   std::map< Simplex_key, int >               idx_to_colorder_; //idx_to_colorder_[simp_key] = order of col in mat
//   std::vector< Simplex_key >                 col_order_; //col_order_[ idx_to_colorder[key] ] = key
//   //------
//   std::map< Simplex_key, Simplex_key >       pairing_; // g -> h, f -> 0 and h -> -1
//   std::map< Simplex_key, Simplex_key >       low_to_col_; // idx -> index of reduced col with // lowest 1 at index idx
//   std::map< Simplex_key, Simplex_key >       col_to_low_; //inverse
//   std::map< Simplex_key, int >               colidx_to_birth_; //return the birth of a column from its idx
// //  std::map< Simplex_key, int >               colidx_to_death_;// same with death






// /*
// * check that everything is well updated ->colidx_to_birth, low_to_col etc  ?
// *
// */
// void forward_arrow( Simplex_handle zzsh )
// {
//   birth_vector_.add_birth_forward();

//   //store the boundary in a set
//   std::set< Simplex_key > col_bsh;
//   for( auto b_sh : cpx_->boundary_simplex_range(zzsh) ) 
//   {  col_bsh.insert( cpx_->key(b_sh) );  } //compute col for boundary of sigma, repr by a set

//   Simplex_key low_idx     = *(col_bsh.rbegin());   //idx of lowest element in col_bsh
//   Simplex_key col_low     =  low_to_col_[low_idx]; //idx of reduced col with low_idx as lowest index
//   Simplex_key paired_idx  =  pairing_[col_low];    //col with which col_low is paired
//   std::vector< Simplex_key > chains_in_H;          //for corresponding indices in H
//   std::vector< Simplex_key > chains_in_G;
//   std::pair< typename std::set< Simplex_key >::iterator, bool > res_insert;

//   while( true ) //reduce col_bsh with cycles of G
//   {
//     if( paired_idx > 0 ) //i.e. col_low \in G and paired_idx \in H
//     { 
//       chains_in_H.push_back(paired_idx); //remember the h with which g is paired
//       chains_in_G.push_back(col_low);    //remember the g
//       //add the column
//       for(auto &cell : mat_[col_low]->col_) {
//           res_insert = col_bsh.insert(cell.key());
//           if( !res_insert.second ) { col_bsh.erase(res_insert.first); }
//       }

//       if(col_bsh.empty()) //col_bsh==0: [\partial sigma]=0 
//       {                    //sigma creates a cycle made with the col_h, h \in chains_in_H
//         for( auto idx_h : chains_in_H ) { //produce the sum of all col in chains_in_H
//           for(auto &cell : mat_[idx_h]->col_ ) {
//             res_insert = col_bsh.insert(cell.key());
//             if( !res_insert.second ) { col_bsh.erase(res_insert.first); }
//           }
//         }
//         //Create a new column with the new cycle value
//         Column * new_col = new Column();
//         for( auto idx : col_bsh )  //add all indices in col_new with canonical order
//         {  
//           Cell * new_cell = new Cell(idx);
//           // todo      transverse_row
//           new_col->col_.push_back( *new_cell );  
//         }

//         Cell * new_cell = new Cell(cpx_->key(zzsh));
//         // todo      transverse_row

//         new_col->col_.push_back( *new_cell ); // and add sigma <- must have biggest idx of all

//         pairing_[cpx_->key(zzsh)] = 0;      //add index of sigma in F
//         mat_[cpx_->key(zzsh)]    = new_col; //add column to the matrix
//         colidx_to_birth_[cpx_->key(zzsh)] = cpx_->key(zzsh); //IWDP birth == key
//         col_to_low_[cpx_->key(zzsh)]      = cpx_->key(zzsh); //death == key
//         low_to_col_[cpx_->key(zzsh)]      = cpx_->key(zzsh);

//       //add row with cpx_->key(zzsh)
//       //......
//       /*******************************************************************/
//       /*******************************************************************/
//       /**************************** todo ********************************/
//       /*******************************************************************/
//       /*******************************************************************/

//         return;
//       }

//     //col_bsh != 0:
//       low_idx     = *(col_bsh.rbegin());   //idx of lowest element
//       col_low     =  low_to_col_[low_idx]; //reduced col with low_idx as lowest index
//       paired_idx  =  pairing_[col_low];    //col with which col_low is paired
//     }
//     else { //no more columns in G to reduce col_bsh, continue with F \cup G     [\partial sigma] != 0

//         std::vector< Simplex_key > chains_in_F;

//         while(true) {
//           if(paired_idx == 0) { chains_in_F.push_back(col_low); } //col_low \in F
//           else                { chains_in_H.push_back(paired_idx); } //paired_idx \in H

//           //add the col
//           for(auto &cell : mat_[col_low]->col_) {
//             res_insert = col_bsh.insert(cell.key());
//             if( !res_insert.second ) { col_bsh.erase(res_insert.first); }
//           }

//           if(col_bsh.empty()) { break; }

//           low_idx     = *(col_bsh.rbegin());         //idx of lowest element
//           col_low     =  low_to_col_[low_idx];      //reduced col with low_idx as lowest index
//           paired_idx  =  pairing_[col_low];         //col with which col_low is paired
//         }

//     //Surjective weak diamond principle
//         // sort the fs by deaths (the closer from WD the bigger)
//         // the death is equal to the lowest index in the column
//         sort(chains_in_F.begin(),chains_in_F.end(), 
//             [this](Simplex_key k1, Simplex_key k2) {return col_to_low_[k1] < col_to_low_[k2];} );

//         //New value for col_fp: col_fp <- col_f1+...+col_fp \in G
//         Simplex_key fp = *(chains_in_F.rbegin()); //the one that gets cut
//         Column *col_fp = mat_[fp]; //corresponding column
//         chains_in_F.pop_back();    //remove fp from chains_in_F
//         for(auto c_idx : chains_in_F) {  plus_equal_column(col_fp,mat_[c_idx]); }

//         //New column col_m+1 that turns col_fp into a boundary
//         //col_m+1 <- sigma + col_h1 +...+ col_hi  with hi in chains_in_H. Here col_bsh is empty
//         for(auto c_idx : chains_in_H) {  
//             for(auto &cell : mat_[c_idx]->col_) {
//             res_insert = col_bsh.insert(cell.key());
//             if( !res_insert.second ) { col_bsh.erase(res_insert.first); }
//           }
//         }
//         Column * new_col = new Column();
//         for( auto idx : col_bsh )  //add all indices in col_new 
//         {  
//           Cell * new_cell = new Cell(idx);
//           // todo      transverse_row
//           new_col->col_.push_back( *new_cell );  
//         }
//         Cell * new_cell = new Cell(cpx_->key(zzsh));
//         new_col->col_.push_back( *new_cell ); // and add sigma <- must have biggest idx of all
//                                               // todo new row
//         mat_[cpx_->key(zzsh)]     = new_col; //insert col for sigma in matrix
//         pairing_[cpx_->key(zzsh)] = -1;      //which belongs to H
//         pairing_[fp]              = cpx_->key(zzsh); //and is paired with fp now in G
//         col_to_low_[cpx_->key(zzsh)] = cpx_->key(zzsh);
//         low_to_col_[cpx_->key(zzsh)] = cpx_->key(zzsh);

//         //re-pairing of the fi
//         auto cmp_birth_vector = [this](Simplex_key k1, Simplex_key k2)->bool 
//                                     {return birth_vector_[k1] > birth_vector_[k2];};//ordered by reverse <=_b
//         std::map< Simplex_key, Simplex_key, decltype(cmp_birth_vector) > birth_to_idx(cmp_birth_vector); 

//         //for f1 to f_{p-1} sorted by <=_d
//         for(auto f_idx : chains_in_F) { birth_to_idx[ colidx_to_birth_[f_idx] ] = f_idx;} //inverse map
//         birth_to_idx[ colidx_to_birth_[fp] ] = fp; //contains p elements, the birth of fp must be available to others

//         auto bti_it = birth_to_idx.begin();  //points to max birth
//         Simplex_key max_birth = birth_to_idx.begin()->first;

//         Simplex_key curr_birth;
//         for(auto f_idx : chains_in_F) {
//           //colidx_to_birth_[f_idx] is the original birth of col_fidx
//           //find which reduced col has this birth
//           bti_it = birth_to_idx.find(colidx_to_birth_[f_idx]); 

//           if(bti_it->first == max_birth) { ++bti_it; } //if its max_birth, give next maximal birth
//           //while the col with curr_birth appears before in sorted chains_in_F
//           while( col_to_low_[bti_it->second] < col_to_low_[f_idx] ) 
//           {
//             plus_equal_column(mat_[f_idx],mat_[bti_it->second]);
//             ++bti_it;
//           }
//           bti_it->second = f_idx;
//         }
//     //update colidx_to_birth
//         bti_it = birth_to_idx.begin(); ++bti_it;
//         for(; bti_it != birth_to_idx.end(); ++bti_it) {
//           colidx_to_birth_[bti_it->second] = bti_it->first;
//         }
//         colidx_to_birth_.erase(fp); //fp not in F anymore

//     //update persistence diagram
//         persistence_diagram.emplace_back(max_birth,cpx_->key(zzsh)); //open interval
//         return;
//     }
//   }  
// }





    // /* Swap the columns .col_ stored in two matrix_chains. */
    // void col_swap(matrix_chain * curr_col, matrix_chain * other_col)
    // {
    //   std::swap(curr_col->column() , other_col->column() );
    //   for(auto &cell : curr_col->column()->col_)  { cell.self_chain_ = curr_col; }
    //   for(auto &cell : other_col->column()->col_) { cell.self_chain_ = other_col;}
    // }




#ifdef _VERBATIM_ZIGZAG_PERSISTENCE_ 
  std::cout << "+++Explicit boundary of inserted simplex: ";
  for(auto v : cpx_->simplex_vertex_range(zzsh)) { std::cout << v << " "; }
  std::cout << std::endl;
  for( auto b_sh : cpx_->boundary_simplex_range(zzsh) )
  {
    std::cout << "  - ";
    std::cout << "key= " << cpx_->key(b_sh) << "  | ";    
    for(auto v : cpx_->simplex_vertex_range(b_sh)) { std::cout << v << " "; }
    std::cout << "\n";
  }
  std::cout << " ===== end boundary\n";
  std::cout << "   boundary: ";
  for(auto idx : col_bsh) {std::cout << idx << " ";}
  std::cout << "\n";
#endif




#ifdef _VERBATIM_ZIGZAG_PERSISTENCE_
  std::cout << "   size boundary = " << col_bsh.size() << "\n";
  std::cout << "   all of lowidx_to_matidx_:\n";
  for(auto pp : lowidx_to_matidx_) { std::cout << pp.first << " "; }
    std::cout << "\n";
#endif



#ifdef _VERBATIM_ZIGZAG_PERSISTENCE_
  std::cout << "** Columns in F: none ";
  std::cout << "\n";
  std::cout << "** Columns in G: ";
  for(auto ch : chains_in_G) { std::cout << ch->lowest_idx_ << " ; ";}
  std::cout << "\n";
  std::cout << "** Columns in H: ";
  for(auto ch : chains_in_H) { std::cout << ch->lowest_idx_ << " ; ";}
  std::cout << "\n";
#endif


  // std::cout << "Columns in G: \n";
  // for(auto col : chains_in_G)
  // {
  //   std::cout << " : [k=" << col->lowest_idx_ <<"] ";
  //   std::cout << "[b=" << col->birth_ <<"] ";
  //   if(col->paired_col_ != nullptr) { std::cout << "[pc=" << col->paired_col_->lowest_idx_ <<"]     ";}
  //   else { std::cout << "[pc=" << "F" <<"]     ";}
  //   std::cout <<"| ";
  //   for(auto &cell : *(col->column())) { std::cout << cell.key() << " "; }
  //   std::cout << " |" << std::endl;
  // }

      
#ifdef _VERBATIM_ZIGZAG_PERSISTENCE_ 
  std::cout << "** Columns in F: ";
  for(auto ch : chains_in_F) { std::cout << ch->lowest_idx_ << " ; ";}
  std::cout << "\n";
  std::cout << "** Columns in G: ";
  for(auto ch : chains_in_G) { std::cout << ch->lowest_idx_ << " ; ";}
  std::cout << "\n";
  std::cout << "** Columns in H: ";
  for(auto ch : chains_in_H) { std::cout << ch->lowest_idx_ << " ; ";}
  std::cout << "\n";
#endif
