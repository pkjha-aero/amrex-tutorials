
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {
        amrex::Print() << "Hello world from AMReX version " << amrex::Version() << "\n";
    }
    {
		Print() << std::endl << ">>>>>> ParallelDescriptor: >>>>>>>>>>>>" << std::endl;
		int myproc = ParallelDescriptor::MyProc(); // Return the rank
		int nprocs = ParallelDescriptor::NProcs(); // Return the number of processes
		if (ParallelDescriptor::IOProcessor()) {
		// Only the I/O process executes this
		}
		int ioproc = ParallelDescriptor::IOProcessorNumber(); // I/O rank
		Print() << "myproc: " << myproc << ", nprocs: "  << nprocs << ", ioproc: "  << ioproc << std::endl;
		ParallelDescriptor::Barrier();
		Print() << "Broadcast 100 ints from the I/O Processor " << std::endl;
		Vector<int> a(100);
		ParallelDescriptor::Bcast(a.data(), a.size(), ParallelDescriptor::IOProcessorNumber());
		// See AMReX_ParallelDescriptor.H for many other Reduce functions
		//ParallelDescriptor::ReduceRealSum(x);
	}
	{
		Print() << std::endl << ">>>>>> ParallelContext: >>>>>>>>>>>>" << std::endl;
		//MPI_Comm subCommA = ....;
		//MPI_Comm subCommB = ....;
		// Add a communicator to ParallelContext.
		// After these pushes, subCommB becomes the
		//     "local" communicator.
		//ParallelContext::push(subCommA);
		//ParallelContext::push(subCommB);

		// Get Global and Local communicator (subCommB).
		MPI_Comm globalComm = ParallelContext::CommunicatorAll();
		MPI_Comm localComm  = ParallelContext::CommunicatorSub();

		// Get local number of ranks and global IO Processor Number.
		int localRanks = ParallelContext::NProcsSub();
		int globalIO     = ParallelContext::IOProcessorNumberAll();
		Print() << "localRanks: " << localRanks << ", globalIO: "  << globalIO << std::endl;

		if (ParallelContext::IOProcessorSub()) {
			// Only the local I/O process executes this
		}

		// Translation of global rank to local communicator rank.
		// Returns MPI_UNDEFINED if comms do not overlap.
		//int localRank = ParallelContext::global_to_local_rank(globalrank);

		// Translations of MPI rank IDs using integer arrays.
		// Returns MPI_UNDEFINED if comms do not overlap.
		//ParallelContext::global_to_local_rank(local_array, global_array, n);
		//ParallelContext::local_to_global_rank(global_array, local_array, n);

		// Remove the last added subcommunicator.
		// This would make "subCommA" the new local communicator.
		// Note: The user still needs to free "subCommB".
		//ParallelContext::pop();
	}
	{
		Print() << std::endl << ">>>>>> Print: >>>>>>>>>>>>" << std::endl;

		Real pi = std::atan(1.0)*4.0;
		// Print on rank 3 with precision of 17 digits
		// SetPrecision does not modify cout's floating-point decimal precision setting.
		Print(3).SetPrecision(17) << pi << "\n";

		int oldprec = std::cout.precision(10);
		Print() << pi << "\n";  // Print with 10 digits

		AllPrint() << "Every process prints\n";  // Print on every process

		std::ofstream ofs("myFile.txt", std::ofstream::out);
		Print(ofs) << "Print to a file" << std::endl;
		ofs.close();

		AllPrintToFile("file.") << "Each process appends to its own file (e.g., file.3)\n";
	}
	{
		Print() << std::endl << ">>>>>> ParmParse: >>>>>>>>>>>>" << std::endl;
		ParmParse pp;

		int nsteps = 0;
		pp.query("nsteps", nsteps);
		Print() << "nsteps: " << nsteps << std::endl;  // 1000

		Real dt;
		pp.get("dt", dt);  // runtime error if dt is not in inputs
		Print() << "dt: " << dt << std::endl;  // 0.03

		Vector<int> numcells;
		// The variable name 'numcells' can be different from parameter name 'ncells'.
		pp.getarr("ncells", numcells);
		Print() << "numcells has size: " << numcells.size() << std::endl;  // 3
		Print() << "numcells has these elements: ";
		for (auto &x:numcells)
			Print() << x << "  ";
		Print() << std::endl;

		Vector<Real> xr {-1.0, 1.0};
		if (!pp.queryarr("xrange", xr)) {
			amrex::Print() << "Cannot find xrange in inputs, "
						   << "so the default {-1.0,1.0} will be used\n";
		}
		Print() << "xr has size: " << xr.size() << std::endl;  // 2
		Print() << "xr has these elements: ";
		for (auto &x:xr)
			Print() << x << "  ";
		Print() << std::endl;

		std::string title;
		pp.query("title", title);  // query string
		Print() << "title: " << title << std::endl; 

		ParmParse pph("hydro");  // with prefix 'hydro'
		Real cfl;
		pph.get("cfl", cfl);    // get parameter with prefix
		Print() << "cfl: " << cfl << std::endl; 
	}
    {
		Print() << std::endl << ">>>>>> IntVect (Arithmetic): >>>>>>>>>>>>" << std::endl;
		IntVect iv(AMREX_D_DECL(19, 0, 5));
		Print() << "iv: " << iv << std::endl;
		for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
	 	   amrex::Print() << "iv[" << idim << "] = " << iv[idim] << std::endl;
		}

		//IntVect iv(AMREX_D_DECL(19, 0, 5));
		IntVect iv2(AMREX_D_DECL(4, 8, 0));
		Print() << "iv2: " << iv2 << std::endl; 
		iv += iv2;  // iv is now (23,8,5)
		Print() << "iv+=iv2: " << iv << std::endl;
		iv *= 2;    // iv is now (46,16,10);
		Print() << "iv*=2: " << iv << std::endl;
    }
    {
		Print() << std::endl << ">>>>>> IntVect (Coarsening): >>>>>>>>>>>>" << std::endl;
		IntVect iv(AMREX_D_DECL(127,127,127));
		Print() << "iv: " << iv << std::endl;
		
		IntVect coarsening_ratio(AMREX_D_DECL(2,2,4));
		Print() << "coarsening_ratio: " << coarsening_ratio << std::endl;
		
		iv.coarsen(2);                 // Coarsen each component by 2
		Print() << "iv (coarsened by 2): " << iv << std::endl;
		
		iv.coarsen(coarsening_ratio);  // Component-wise coarsening
		Print() << "iv (coarsened by 'coarsening_ratio'): " << iv << std::endl;
		
		const auto& iv2 = amrex::coarsen(iv, 2); // Return an IntVect w/o modifying iv
		Print() << "iv2 (iv coarsened by 2, but w/o modifying iv): " << iv2 << std::endl;
		
		IntVect iv3 = amrex::coarsen(iv, coarsening_ratio); // iv not modified
		Print() << "iv3 (iv coarsened by 'coarsening_ratio', but w/o modifying iv): " << iv3 << std::endl;
    }
    {
		Print() << std::endl << ">>>>>> Box (Cell/Node-Based): >>>>>>>>>>>>" << std::endl;
		IntVect lo(AMREX_D_DECL(64,64,64));
		IntVect hi(AMREX_D_DECL(127,127,127));
		IndexType typ({AMREX_D_DECL(1,1,1)});
		Print() << "lo: " << lo << std::endl;
		Print() << "hi: " << hi << std::endl;
		Print() << "typ: " << typ << std::endl;
		Box cc(lo,hi);        // By default, Box is cell based.
		Box nd(lo,hi+1,typ);  // Construct a nodal Box.
		Print() << "A cell-centered Box: " << cc << "\n";
		Print() << "An all nodal Box:    " << nd << "\n";
    }
    {
		Print() << std::endl << ">>>>>> Box (Cell/Node Conversion): >>>>>>>>>>>>" << std::endl;
		Box b0 ({64,64,64}, {127,127,127}); // Index type: (cell, cell, cell)
		Box b1 = surroundingNodes(b0);  // A new Box with type (node, node, node)
		
		Print() << "b0: " << b0 << std::endl;                  // Still ((64,64,64) (127,127,127) (0,0,0))
		Print() << "b1 (= surrounding nodes of b0): " << b1 << std::endl;                  // ((64,64,64) (128,128,128) (1,1,1))


		Box b2 = enclosedCells(b1);     // A new Box with type (cell, cell, cell)
		Print() << "b2 (Enclosed cells of b1): " << b2 << std::endl;
		if (b2 == b0) {                 // Yes, they are identical.
		   Print() << "b0 and b2 are identical!\n";
		}

		Box b3 = convert(b0, {0,1,0});  // A new Box with type (cell, node, cell)
		Print() << "b3 (after converting b0 by {0, 1, 0}): " << b3 << std::endl;                  // ((64,64,64) (127,128,127) (0,1,0))

		b3.convert({0,0,1});            // Convert b0 to type (cell, cell, node)
		Print() << "b3 (after converting it by {0, 0, 1}): " << b3 << std::endl;                  // ((64,64,64) (127,127,128) (0,0,1))

		Print() << "Surrounding nodes of b3:" << b3.surroundingNodes() << std::endl;          //  Exercise for you
		Print() << "Enclosed Cells of b3: " << b3.enclosedCells() << std::endl;             //  Exercise for you
    }  
    {
		Print() << std::endl << ">>>>>> Box (Refining and Coarsening): >>>>>>>>>>>>" << std::endl;
		Print() << "Cell-centered: " << std::endl;
		Box ccbx ({16,16,16}, {31,31,31});
		Print() << "ccbx: " << ccbx << std::endl;
		ccbx.refine(2);
		Print() << "ccbx (refined by 2): " << ccbx << std::endl;                   // ((32,32,32) (63,63,63) (0,0,0))
		Print() << "ccbx (coarsened by 2): " << ccbx.coarsen(2) << std::endl;        // ((16,16,16) (31,31,31) (0,0,0))

		Print() << "Node-centered: " << std::endl;
		Box ndbx ({16,16,16}, {32,32,32}, {1,1,1});
		Print() << "ndbx: " << ndbx << std::endl;
		ndbx.refine(2);
		Print() << "ndbx (refined by 2): " << ndbx << std::endl;                   // ((32,32,32) (64,64,64) (1,1,1))
		Print() << "ndbx (coarsened by 2): "<< ndbx.coarsen(2) << std::endl;        // ((16,16,16) (32,32,32) (1,1,1))

		Print() << "Face-centered in x-dir: \n";
		Box facebx ({16,16,16}, {32,31,31}, {1,0,0});
		Print() << "facebx: " << facebx << std::endl;
		facebx.refine(2);
		Print() << "facebx (refined by 2): " << facebx << std::endl; // ((32,32,32) (64,63,63) (1,0,0))
		Print() << "facebx (coarsened by 2): "<< facebx.coarsen(2) << std::endl;      // ((16,16,16) (32,31,31) (1,0,0))
		
		Print() << "Uncoarsenable: \n";
		Box uncoarsenable ({16,16,16}, {30,30,30});
		Print() << uncoarsenable.coarsen(2) << std::endl; // ((8,8,8), (15,15,15));
		Print() << uncoarsenable.refine(2) << std::endl;  // ((16,16,16), (31,31,31));
	                                      // Different from the original!
    }
    {
		Print() << std::endl << ">>>>>> Box (Intersection): >>>>>>>>>>>>" << std::endl;
		Box b0 ({16,16,16}, {31,31,31});
		Box b1 ({ 0, 0,30}, {23,23,63});
		Print() << "b0: " << b0 << std::endl;
		Print() << "b1: " << b1 << std::endl;

		if (b0.intersects(b1)) {                  // true
		    Print() << "b0 and b1 intersect.\n";
		}

		Box b2 = b0 & b1;     // b0 and b1 unchanged
		Print() << "b2 = b0 & b1: " << b2 << std::endl;        // ((16,16,30) (23,23,31) (0,0,0))

		Box b3 = surroundingNodes(b0) & surroundingNodes(b1); // b0 and b1 unchanged
		Print() << "b3 = surroundingNodes(b0) & surroundingNodes(b1): " << b3 << std::endl;        // ((16,16,30) (24,24,32) (1,1,1))

		b0 &= b2;             // b2 unchanged
		Print() << "b0 &= b2: " << b0 << std::endl;        // ((16,16,30) (23,23,31) (0,0,0))
	        
	    Print() << "Attempting intersection of cell-centered and node-centered Box \n";
		//b0 &= b3;             // Runtime error because of type mismatch!	
    }
    {
		Print() << std::endl << ">>>>>> Dim3 and XDim3: >>>>>>>>>>>>" << std::endl;
		IntVect iv(AMREX_D_DECL(64,64,64));	
		Dim3 d3 = iv.dim3();
	    Print() << "d3 (Dim3): " << d3 << std::endl;

		Box bx ({ 0, 0,30}, {23,23,63});
		Print() << "bx: "<< bx << std::endl;	
		Dim3 lo = lbound(bx);
		Dim3 hi = ubound(bx);
		Print() << "lo of bx: " << lo << std::endl;
		Print() << "hi of bx: " << hi << std::endl;
		Print() << "length of bx: " << length(bx);
		Print() << std::endl << std::endl;
    }
    {
		Print() << std::endl << ">>>>>> RealBox and Geometry: >>>>>>>>>>>>" << std::endl;
		int n_cell = 64;

		// This defines a Box with n_cell cells in each direction.
		Box domain(IntVect{AMREX_D_DECL(       0,        0,        0)},
				   IntVect{AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1)});
		Print() << "domain (index): " << domain << std::endl;

		// This defines the physical box, [-1,1] in each direction.
		RealBox real_box({AMREX_D_DECL(-1.0,-1.0,-1.0)},
						 {AMREX_D_DECL( 1.0, 1.0, 1.0)});
		Print() << "real_box (real coordinates): " << real_box << std::endl;
		
		// This says we are using Cartesian coordinates
		int coord = 0;

		// This sets the boundary conditions to be doubly or triply periodic
		Array<int,AMREX_SPACEDIM> is_periodic {AMREX_D_DECL(1,1,1)};
		Print() << "is_periodic: " << is_periodic << std::endl;

		// This defines a Geometry object
		Geometry geom(domain, real_box, coord, is_periodic);
		Print() << "geom: " << geom << std::endl;

		const auto problo = geom.ProbLoArray(); // Lower corner of the physical
												// domain.  The return type is
												// GpuArray<Real,AMREX_SPACEDIM>.
		Print() << "problo: " << *geom.ProbLo() << std::endl;
		
		Real yhi = geom.ProbHi(1);              // y-direction upper corner
		Print() << "yhi: " << yhi << std::endl;
		
		const auto dx = geom.CellSizeArray();   // Cell size for each direction.
		//Print() << "dx: " << dx << std::endl;
		
		const Box& domain_ret = geom.Domain();      // Index domain
		Print() << "domain (retrieved from geometry): " << domain_ret << std::endl;
		
		bool is_per = geom.isPeriodic(0);       // Is periodic in x-direction?
		Print() << "is_per: " << geom.isPeriodic() << std::endl;
		Print() << "is_per (0): " << is_per << std::endl;
		
		//if (geom.isAllPeriodic()) {}            // Periodic in all direction?
		//if (geom.isAnyPeriodic()) {}            // Periodic in any direction?
	}
    Print() << std::endl << ">>>>>> GOOD BYE >>>>>>>>>>>>" << std::endl;
    amrex::Finalize();
}
