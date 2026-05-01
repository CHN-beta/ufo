{
  inputs.nixos.url = "github:CHN-beta/nixos";

  outputs = inputs:
  {
    packages.x86_64-linux = let inherit (inputs.nixos.packages.x86_64-linux) pkgs; in rec
    {
      default = pkgs.pkgsStatic.localPkgs.buildUfo ./.;
      release = pkgs.runCommand "x86_64-linux-musl-static.zip" {}
      ''
        ${pkgs.zip}/bin/zip -r $out ${default}/bin/ufo
      '';
    };
  };
}
