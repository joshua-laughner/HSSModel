using Documenter, HSSModel;

makedocs(sitename="HSSModel");
if length(ARGS) > 0 && ARGS[1] == "deploy"
    deploydocs(repo = "github.com/joshua-laughner/HSSModel.git",)
else
    println("To deploy to gh-pages, add `deploy` as a command line argument")
end
