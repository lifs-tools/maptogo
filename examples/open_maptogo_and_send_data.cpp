// open_maptogo_and_send_data.cpp
// Sends large data to a Flask/Dash endpoint via HTTP POST
//
// Dependencies (header-only, no build system needed):
//   cpp-httplib:      https://github.com/yhirose/cpp-httplib  (single header: httplib.h)
//   nlohmann/json:    https://github.com/nlohmann/json        (single header: json.hpp)
//
// Compile:
//   g++ -std=c++17 -O3 open_maptogo_and_send_data.cpp -o open_maptogo_and_send_data
//   # On Linux you may also need: -lpthread
//
// Usage:
//   ./open_maptogo_and_send_data

#include "httplib.h"   // cpp-httplib (drop next to this file)
#include "json.hpp"    // nlohmann/json (drop next to this file)

#include <iostream>
#include <vector>
#include <string>

using json = nlohmann::json;


// ---------------------------------------------------------------------------
// Send to Flask
// ---------------------------------------------------------------------------
std::string post_data(const std::string& host, int port,
               const std::string& endpoint,
               const std::string& body) {

    httplib::Client cli(host, port);
    cli.set_connection_timeout(5);
    cli.set_read_timeout(30);
    cli.set_write_timeout(30);

    auto res = cli.Post(endpoint.c_str(), body, "application/json");

    if (!res) {
        return "";
    }

    if (res->status == 200) {
        return res->body;
    } else {
        return "";
    }
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main() {
    const std::string HOST     = "localhost";
    const int         PORT     = 8040;
    const std::string ENDPOINT = "/submit";

    json data;
    // if separate_updown_switch is true, [up/down]regulated_[l/p/m/t] fields with be
    // taken, else only regulated_[l/p/m/t] fields
    data["separate_updown_switch"] = "true";
    data["all_lipids"] = {"PS 34:2", "PS 34:1", "PS 36:4", "PS 36:3", "PS 36:2", "PS 36:1", "PS 38:5", "PS 38:4", "PS 38:3", "PS 38:2", "PS 38:1", "PS 39:3", "PS 40:7", "PS 40:6", "PS 40:5", "PS 40:4", "PS 40:3", "PS 40:2", "PS 40:1", "PI 32:1", "PI 34:2", "PI 34:1", "PI 36:5", "PI 36:4", "PI 36:3", "PI 36:2", "PI 37:4", "PI 37:3", "PI 37:2", "PI 38:6", "PI 38:5", "PI 38:4", "PI 38:3", "PI 39:5", "PI 40:7", "PI 40:6", "PI 40:5", "TAG 44:1", "TAG 46:2", "TAG 47:2", "TAG 48:1", "TAG 48:2", "TAG 48:3", "TAG 49:2", "TAG 49:3", "TAG 50:2", "TAG 50:3", "TAG 50:4", "TAG 51:2", "TAG 51:3", "TAG 52:3", "TAG 52:4", "TAG 52:5", "TAG 53:3", "TAG 53:4", "TAG 54:3", "TAG 54:4", "TAG 54:5", "TAG 54:6", "TAG 54:7", "TAG 56:2", "TAG 56:3", "TAG 56:4", "TAG 56:5", "TAG 56:6", "TAG 56:7", "TAG 56:8", "TAG 58:2", "TAG 58:3", "TAG 58:4", "TAG 58:5", "TAG 58:6", "TAG 58:7", "TAG 58:8", "TAG 58:9", "TAG 60:3", "TAG 60:4", "PC-O 32:0", "PC-O 32:1", "PC-O 34:0", "PC-O 34:1", "PC-O 36:3", "PC-O 38:3", "PC-O 38:4", "PC 32:0", "PC 32:1", "PC 32:2", "PC 33:1", "PC 33:2", "PC 34:1", "PC 34:2", "PC 34:3", "PC 35:1", "PC 35:2", "PC 36:1", "PC 36:2", "PC 36:3", "PC 36:4", "PC 38:4", "PC 38:5", "PC 40:5", "PC 40:6", "LPC 14:0", "LPC 15:0", "LPC 16:1", "LPC 16:0", "LPC 17:0", "LPC 18:2", "LPC 18:1", "LPC 18:0", "LPC 20:4", "LPC 20:3", "LPC 20:1", "LPC 20:0", "LPC 22:6", "LPC 22:5", "LPC 22:4", "DG 30:0", "DG 32:0", "DG 32:1", "DG 34:0", "DG 34:1", "DG 34:2", "DG 36:1", "DG 36:2", "DG 36:4", "DG 38:5", "DG 43:6", "DG 46:1", "DG 48:1", "DG 50:1", "Cer 34:2", "Cer 34:1", "Cer 36:1", "Cer 40:2", "Cer 40:1", "Cer 42:3", "Cer 42:2", "Cer 42:1", "SM 32:1", "SM 33:1", "SM 34:1", "SM 34:0", "SM 35:1", "SM 36:1", "SM 38:2", "SM 38:1", "SM 39:1", "SM 40:2", "SM 40:1", "SM 41:1", "SM 42:2", "SM 42:1", "PG 32:0", "PG 32:1", "PG 34:1", "PG 34:2", "PG 34:3", "PG 34:4", "PG 36:1", "PG 36:3", "PG 36:4", "PG 38:2", "PG 38:3", "PG 38:4", "PG 38:5", "PG 38:6", "PG 40:3", "PG 40:5", "PG 40:6", "PG 40:7", "PG 42:7", "PE-O 32:2", "PE-O 34:1", "PE-O 34:2", "PE-O 34:3", "PE-O 35:2", "PE-O 35:4", "PE-O 36:1", "PE-O 36:2", "PE-O 36:3", "PE-O 36:4", "PE-O 36:5", "PE-O 36:6", "PE-O 37:4", "PE-O 37:5", "PE-O 37:6", "PE-O 38:1", "PE-O 38:2", "PE-O 38:3", "PE-O 38:4", "PE-O 38:5", "PE-O 38:6", "PE-O 38:7", "PE-O 40:2", "PE-O 40:3", "PE-O 40:5", "PE 32:1", "PE 32:2", "PE 33:1", "PE 34:1", "PE 34:2", "PE 35:1", "PE 35:2", "PE 36:1", "PE 36:2", "PE 36:3", "PE 36:4", "PE 36:6", "PE 37:2", "PE 37:3", "PE 37:4", "PE 38:2", "PE 38:3", "PE 38:4", "PE 38:5", "PE 38:6", "PE 38:7", "PE 39:5", "PE 39:6", "PE 40:1", "PE 40:2", "PE 40:4", "PE 40:5", "PE 40:6", "PE 40:7", "PE 42:7", "PA 32:0", "PA 32:1", "PA 34:0", "PA 34:1", "PA 34:2", "PA 35:2", "PA 36:2", "PA 36:1", "PA 38:4", "LPE 16:0", "LPE 16:1", "LPE 17:0", "LPE 18:0", "LPE 18:1", "LPE 20:4", "LPE 20:5", "LPE 22:5", "LPE 22:6", "LPA 16:0", "LPA 18:3", "LPA 18:0"};
    data["upregulated_lipids"] = {"PS 34:2", "PS 34:1", "PS 36:4", "PS 36:3", "PS 36:2", "PS 36:1", "PS 38:5", "PS 38:4", "PS 38:2", "PS 38:1", "PS 39:3", "PS 40:7", "PS 40:6", "PS 40:5", "PS 40:4", "PS 40:3", "PS 40:2", "PS 40:1", "PI 32:1", "PI 34:2", "PI 34:1", "PI 36:5", "PI 36:4", "PI 36:3", "PI 36:2", "PI 37:4", "PI 37:3", "PI 37:2", "PI 38:6", "PI 38:5", "PI 38:4", "PI 38:3", "PI 39:5", "PI 40:7", "PI 40:6", "PI 40:5", "TAG 44:1", "TAG 46:2", "TAG 47:2", "TAG 48:1", "TAG 48:2", "TAG 49:2", "TAG 49:3", "TAG 50:2", "TAG 52:5", "TAG 54:7", "TAG 58:2", "PC 32:1", "PC 32:2", "PC 33:1", "PC 33:2", "PC 34:1", "PC 34:2", "PC 35:1", "PC 36:1", "LPC 14:0", "LPC 15:0", "LPC 16:1", "LPC 16:0", "LPC 17:0", "LPC 18:2", "LPC 18:1", "LPC 18:0", "LPC 20:3", "LPC 22:6", "DG 30:0", "DG 32:0", "DG 32:1", "DG 34:0", "DG 34:1", "DG 46:1", "Cer 34:2", "Cer 40:2", "Cer 42:3", "SM 33:1", "PG 32:0", "PG 32:1", "PG 34:2", "PG 34:3", "PE-O 32:2", "PE-O 36:2", "PE 32:1", "PE 32:2", "PE 33:1", "PE 34:1", "PE 34:2", "PE 36:4", "PE 36:6", "PE 38:7", "PE 40:6", "PA 32:0", "PA 32:1", "PA 34:0", "PA 34:1", "PA 34:2", "PA 35:2", "LPE 16:0", "LPE 17:0", "LPE 18:0", "LPE 18:1", "LPA 16:0", "LPA 18:0"};
    data["downregulated_lipids"] =  {"TAG 58:5", "TAG 58:6", "PC-O 34:1", "PC-O 36:3", "PC-O 38:3", "PC-O 38:4", "DG 36:2", "SM 34:0", "PG 34:4", "PG 36:1", "PG 36:3", "PG 36:4", "PG 38:2", "PG 38:3", "PG 38:4", "PG 38:5", "PG 40:3", "PG 40:5", "PG 40:6", "PG 40:7", "PG 42:7", "PE-O 34:2", "PE-O 35:4", "PE-O 36:3", "PE-O 36:4", "PE-O 36:5", "PE-O 37:5", "PE-O 38:1", "PE-O 38:2", "PE-O 38:3", "PE-O 38:4", "PE-O 38:5", "PE-O 40:3", "PE-O 40:5", "PE 37:2", "PE 38:2", "PE 38:3", "PE 39:5", "PE 40:1", "PE 40:2", "PE 40:4", "PE 40:7", "PE 42:7", "PA 36:1", "PA 38:4", "LPE 20:5", "LPE 22:5", "LPE 22:6", "LPA 18:3"};
    // data["regulated_lipids"] = ...
    // data["all_proteins"] = ...
    // data["regulated_proteins"] = ...
    // data["upregulated_proteins"] = ...
    // data["downregulated_proteins"] = ...
    // data["all_metabolites"] = ...
    // data["regulated_metabolites"] = ...
    // data["upregulated_metabolites"] = ...
    // data["downregulated_metabolites"] = ...
    // data["all_transcripts"] = ...
    // data["regulated_transcripts"] = ...
    // data["upregulated_transcripts"] = ...
    // data["downregulated_transcripts"] = ...
    // data["use_bounded_fatty_acyls"] = ...
    data["use_lipids"] = true;
    // data["use_proteins"] = ...
    // data["use_metabolites"] = ...
    // data["use_transcripts"] = ...
    // data["organism"] = ...
    // data["molecule_handling"] = ...
    // data["regulated_molecule_handling"] = ...
    // data["term_representation"] = ...
    // data["test_method"] = ...
    // data["domains"] = ...
    // data["use_upregulated_lipids"] = ...
    // data["use_downregulated_lipids"] = ...
    // data["use_upregulated_proteins"] = ...
    // data["use_downregulated_proteins"] = ...
    // data["use_upregulated_metabolites"] = ...
    // data["use_downregulated_metabolites"] = ...
    // data["use_upregulated_transcripts"] = ...
    // data["use_downregulated_transcripts"] = ...

    std::string body = data.dump();
    std::string response = post_data(HOST, PORT, ENDPOINT, body);

    if (!response.empty()) {
        auto j = json::parse(response);
        std::string url = "http://localhost:8040";
        if (j.contains("uid")) {
            std::string uid = j["uid"];
            url += "/?uid=" + uid;
        }

        #ifdef __linux__
            int _ = system(("xdg-open \"" + url + "\"").c_str());
        #elif __APPLE__
            int _ = system(("open \"" + url + "\"").c_str());
        #elif _WIN32
            int _ = system(("start \"" + url + "\"").c_str());
        #endif

        return 0;
    }

    return 1;
}
