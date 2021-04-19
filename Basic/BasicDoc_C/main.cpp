
#include <AMReX.H>
#include <AMReX_Print.H>
//#include <AMReX_IntVect.H>
#include <AMReX_Box.H>

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {
        amrex::Print() << "Hello world from AMReX version " << amrex::Version() << "\n";
    }
    {
		IntVect iv(AMREX_D_DECL(19, 0, 5));
		amrex::Print() << iv << std::endl;
		for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
	 	   amrex::Print() << "iv[" << idim << "] = " << iv[idim] << "\n";
		}

		//IntVect iv(AMREX_D_DECL(19, 0, 5));
		IntVect iv2(AMREX_D_DECL(4, 8, 0));
		iv += iv2;  // iv is now (23,8,5)
		Print() << iv << std::endl;
		iv *= 2;    // iv is now (46,16,10);
		Print() << iv << std::endl;
    }
    {
		IntVect iv(AMREX_D_DECL(127,127,127));
		Print() << iv << std::endl;
		IntVect coarsening_ratio(AMREX_D_DECL(2,2,4));
		iv.coarsen(2);                 // Coarsen each component by 2
		Print() << iv << std::endl;
		iv.coarsen(coarsening_ratio);  // Component-wise coarsening
		Print() << iv << std::endl;
		const auto& iv2 = amrex::coarsen(iv, 2); // Return an IntVect w/o modifying iv
		Print() << iv2 << std::endl;
		IntVect iv3 = amrex::coarsen(iv, coarsening_ratio); // iv not modified
		Print() << iv3 << std::endl;
    }
    // BOX
    {
		IntVect lo(AMREX_D_DECL(64,64,64));
		IntVect hi(AMREX_D_DECL(127,127,127));
		IndexType typ({AMREX_D_DECL(1,1,1)});
		Box cc(lo,hi);        // By default, Box is cell based.
		Box nd(lo,hi+1,typ);  // Construct a nodal Box.
		Print() << "A cell-centered Box " << cc << "\n";
		Print() << "An all nodal Box    " << nd << "\n";
    }
    {
		Box b0 ({64,64,64}, {127,127,127}); // Index type: (cell, cell, cell)

		Box b1 = surroundingNodes(b0);  // A new Box with type (node, node, node)
		Print() << b1 << std::endl;                  // ((64,64,64) (128,128,128) (1,1,1))
		Print() << b0 << std::endl;                  // Still ((64,64,64) (127,127,127) (0,0,0))

		Box b2 = enclosedCells(b1);     // A new Box with type (cell, cell, cell)
		if (b2 == b0) {                 // Yes, they are identical.
		   Print() << "b0 and b2 are identical!\n";
		}

		Box b3 = convert(b0, {0,1,0});  // A new Box with type (cell, node, cell)
		Print() << b3 << std::endl;                  // ((64,64,64) (127,128,127) (0,1,0))

		b3.convert({0,0,1});            // Convert b0 to type (cell, cell, node)
		Print() << b3 << std::endl;                  // ((64,64,64) (127,127,128) (0,0,1))

		Print() << "Surrounding nodes:" << b3.surroundingNodes() << std::endl;          //  Exercise for you
		Print() << "Enclosed Cells: " << b3.enclosedCells() << std::endl;             //  Exercise for you
    }
    /*{
	const IntVect& smallEnd () const&;  // Get the small end of the Box
	int bigEnd (int dir) const;         // Get the big end in dir direction
	const int* loVect () const&;        // Get a const pointer to the lower end
	const int* hiVect () const&;        // Get a const pointer to the upper end
    } */  
    {
		Print() << std::endl << "Refining and Coarsening of Box" << std::endl; 
		Box ccbx ({16,16,16}, {31,31,31});
		ccbx.refine(2);
	        Print() << "Cell-centered: \n";
		Print() << ccbx << std::endl;                   // ((32,32,32) (63,63,63) (0,0,0))
		Print() << ccbx.coarsen(2) << std::endl;        // ((16,16,16) (31,31,31) (0,0,0))

		Box ndbx ({16,16,16}, {32,32,32}, {1,1,1});
		ndbx.refine(2);
		Print() << "Node-centered: \n";
		Print() << ndbx << std::endl;                   // ((32,32,32) (64,64,64) (1,1,1))
		Print() << ndbx.coarsen(2) << std::endl;        // ((16,16,16) (32,32,32) (1,1,1))

		Box facebx ({16,16,16}, {32,31,31}, {1,0,0});
		facebx.refine(2);
		Print() << "Face-centered in x-dir: \n";
		Print() << facebx << std::endl;                 // ((32,32,32) (64,63,63) (1,0,0))
		Print() << facebx.coarsen(2) << std::endl;      // ((16,16,16) (32,31,31) (1,0,0))
		
		Print() << "Uncoarsenable: \n";
		Box uncoarsenable ({16,16,16}, {30,30,30});
		Print() << uncoarsenable.coarsen(2) << std::endl; // ((8,8,8), (15,15,15));
		Print() << uncoarsenable.refine(2) << std::endl;  // ((16,16,16), (31,31,31));
	                                      // Different from the original!
    }
    /*                                     
    {                                     //
	bool contains (const Box& b) const;
	bool strictly_contains (const Box& b) const;
	bool contains (const IntVect& p) const;
	bool strictly_contains (const IntVect& p) const;
    }*/
    {
		Print() << std::endl << "Intersection of boxes" << std::endl;
		Box b0 ({16,16,16}, {31,31,31});
		Box b1 ({ 0, 0,30}, {23,23,63});
		if (b0.intersects(b1)) {                  // true
		    Print() << "b0 and b1 intersect.\n";
		}

		Box b2 = b0 & b1;     // b0 and b1 unchanged
		Print() << b2 << std::endl;        // ((16,16,30) (23,23,31) (0,0,0))

		Box b3 = surroundingNodes(b0) & surroundingNodes(b1); // b0 and b1 unchanged
		Print() << b3 << std::endl;        // ((16,16,30) (24,24,32) (1,1,1))

		b0 &= b2;             // b2 unchanged
		Print() << b0 << std::endl;        // ((16,16,30) (23,23,31) (0,0,0))
	        
	        Print() << "Attempting intersection of cell-centered and node-centered Box \n";
		//b0 &= b3;             // Runtime error because of type mismatch!	
    }
    {
		Print() << std::endl << "Dim3 and XDim3" << std::endl;
		IntVect iv(AMREX_D_DECL(64,64,64));	
		Dim3 d3 = iv.dim3();
	        Print() << "Dim3" << d3 << std::endl;

		Box bx ({ 0, 0,30}, {23,23,63});
		Print() << "bx: "<< bx << std::endl;	
		Dim3 lo = lbound(bx);
		Dim3 hi = ubound(bx);
		Print() << "lo of bx: " << lo << std::endl;
		Print() << "hi of bx: " << hi << std::endl;
		Print() << "length of bx: " << length(bx);
		Print() << std::endl << std::endl;
    }
    amrex::Finalize();
}
