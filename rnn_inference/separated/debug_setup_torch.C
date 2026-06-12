// setup_torch.C - macOS compatible version
// Usage: root -l -q setup_torch.C

void setup_torch() {
    std::cout << "=== Setting up LibTorch for ROOT (macOS) ===" << std::endl;
    
    // Get LibTorch path from environment
    TString torch_root = gSystem->Getenv("TORCH_ROOT");
    if (torch_root == "") {
        std::cerr << "TORCH_ROOT environment variable not set!" << std::endl;
        std::cerr << "Please set it first:" << std::endl;
        std::cerr << "export TORCH_ROOT=/path/to/libtorch" << std::endl;
        std::cerr << "export DYLD_LIBRARY_PATH=$TORCH_ROOT/lib:$DYLD_LIBRARY_PATH" << std::endl;
        return;
    }
    
    std::cout << "LibTorch path: " << torch_root << std::endl;
    
    // Check if LibTorch directory exists
    if (!gSystem->AccessPathName(torch_root)) {
        std::cout << "LibTorch directory found" << std::endl;
    } else {
        std::cerr << "LibTorch directory not found: " << torch_root << std::endl;
        return;
    }
    
    // Show what's in the lib directory for debugging
    std::cout << "Contents of lib directory:" << std::endl;
    gSystem->Exec(Form("ls -la %s/lib/", torch_root.Data()));
    
    // Add include paths - torch.h is in csrc/api/include subdirectory
    TString include_path1 = Form("-I%s/include", torch_root.Data());
    TString include_path2 = Form("-I%s/include/torch/csrc/api/include", torch_root.Data());
    
    gSystem->AddIncludePath(include_path1);
    gSystem->AddIncludePath(include_path2);
    std::cout << "✓ Include paths added:" << std::endl;
    std::cout << "  " << include_path1 << std::endl;
    std::cout << "  " << include_path2 << std::endl;
    
    // Add library path
    TString lib_path = Form("-L%s/lib", torch_root.Data());
    gSystem->AddLinkedLibs(lib_path);
    std::cout << "✓ Library path added" << std::endl;
    
    // Try both .dylib and .so extensions
    std::vector<TString> extensions = {".dylib", ".so"};
    TString c10_lib, torch_cpu_lib, torch_lib;
    
    for (const auto& ext : extensions) {
        c10_lib = Form("%s/lib/libc10%s", torch_root.Data(), ext.Data());
        if (!gSystem->AccessPathName(c10_lib)) {
            std::cout << "Found libc10 with extension: " << ext << std::endl;
            break;
        }
    }
    
    for (const auto& ext : extensions) {
        torch_cpu_lib = Form("%s/lib/libtorch_cpu%s", torch_root.Data(), ext.Data());
        if (!gSystem->AccessPathName(torch_cpu_lib)) {
            std::cout << "Found libtorch_cpu with extension: " << ext << std::endl;
            break;
        }
    }
    
    for (const auto& ext : extensions) {
        torch_lib = Form("%s/lib/libtorch%s", torch_root.Data(), ext.Data());
        if (!gSystem->AccessPathName(torch_lib)) {
            std::cout << "Found libtorch with extension: " << ext << std::endl;
            break;
        }
    }
    
    // Check if libraries exist
    if (gSystem->AccessPathName(c10_lib)) {
        std::cerr << "Library not found: " << c10_lib << std::endl;
        std::cerr << "Available files in lib/:" << std::endl;
        gSystem->Exec(Form("ls %s/lib/ | grep -E '(c10|torch)'", torch_root.Data()));
        return;
    }
    
    try {
        std::cout << "Loading libraries in dependency order..." << std::endl;
        
        // Load global dependencies first
        TString global_deps = Form("%s/lib/libtorch_global_deps.so", torch_root.Data());
        if (!gSystem->AccessPathName(global_deps)) {
            if (gSystem->Load(global_deps) < 0) {
                std::cout << "Warning: Could not load global deps, continuing..." << std::endl;
            } else {
                std::cout << "Loaded: libtorch_global_deps.so" << std::endl;
            }
        }
        
        // Try loading c10 with full path and dependencies
        std::cout << "Attempting to load c10..." << std::endl;
        if (gSystem->Load(c10_lib) < 0) {
            std::cerr << "Failed to load c10. Trying alternative approach..." << std::endl;
            
            // Alternative: Use dlopen approach
            gSystem->AddLinkedLibs("-lc10");
            std::cout << "Added c10 to linked libs" << std::endl;
        } else {
            std::cout << "Loaded: " << c10_lib << std::endl;
        }
        
        // Load torch_cpu
        std::cout << "Attempting to load torch_cpu..." << std::endl;
        if (gSystem->Load(torch_cpu_lib) < 0) {
            std::cerr << "Failed to load torch_cpu. Trying alternative..." << std::endl;
            gSystem->AddLinkedLibs("-ltorch_cpu");
            std::cout << "Added torch_cpu to linked libs" << std::endl;
        } else {
            std::cout << "Loaded: " << torch_cpu_lib << std::endl;
        }
        
        // Load main torch library
        std::cout << "Attempting to load torch..." << std::endl;
        if (gSystem->Load(torch_lib) < 0) {
            std::cerr << "Failed to load torch. Trying alternative..." << std::endl;
            gSystem->AddLinkedLibs("-ltorch");
            std::cout << "Added torch to linked libs" << std::endl;
        } else {
            std::cout << "Loaded: " << torch_lib << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Exception loading libraries: " << e.what() << std::endl;
        std::cout << "Trying fallback approach with linked libs..." << std::endl;
        
        // Fallback: add all as linked libraries
        gSystem->AddLinkedLibs("-lc10 -ltorch_cpu -ltorch");
        std::cout << "Added libraries via linking" << std::endl;
    }
    
    std::cout << "LibTorch setup complete!" << std::endl;
    std::cout << std::endl;
    std::cout << "Now you can test:" << std::endl;
    std::cout << "#include <torch/torch.h>" << std::endl;
    std::cout << "auto tensor = torch::ones({2, 3});" << std::endl;
    std::cout << "tensor.print();" << std::endl;
}
