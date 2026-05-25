#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "implot.h"

#include "Tensor.h"
#include "Solver.h"

#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <limits>

using namespace std;

constexpr double PI = 3.14159265358979323846;

static const glm::vec3 OrbitColors[5] = {
    {0.20f, 0.60f, 1.00f},
    {1.00f, 0.45f, 0.20f},
    {0.30f, 0.90f, 0.55f},
    {0.95f, 0.85f, 0.20f},
    {0.85f, 0.30f, 0.80f},
};

glm::vec3 ToCartesian(const StateVector& X, double a)
{
    double r     = X.X[1];
    double theta = X.X[2];
    double phi   = X.X[3];
    double rho   = sqrt(r * r + a * a);
    return glm::vec3(
        (float)(rho * sin(theta) * cos(phi)),
        (float)(rho * sin(theta) * sin(phi)),
        (float)(r   * cos(theta))
    );
}

float  yaw        = -90.f;
float  pitch      =   0.f;
float  radius     =  50.f;
double lastX      = 500.0;
double lastY      = 500.0;
bool   firstMouse = true;
bool   mouseActive= true;

void mouse_callback(GLFWwindow*, double xpos, double ypos)
{
    if (!mouseActive) return;
    if (firstMouse) { lastX = xpos; lastY = ypos; firstMouse = false; }
    float dx = (float)(xpos - lastX);
    float dy = (float)(lastY - ypos);
    lastX = xpos; lastY = ypos;
    yaw   += dx * 0.1f;
    pitch += dy * 0.1f;
    pitch  = glm::clamp(pitch, -89.f, 89.f);
}

void scroll_callback(GLFWwindow*, double, double yoffset)
{
    radius -= (float)yoffset;
    radius  = glm::clamp(radius, 5.f, 300.f);
}

GLuint CompileShader(GLenum type, const char* src)
{
    GLuint s = glCreateShader(type);
    glShaderSource(s, 1, &src, nullptr);
    glCompileShader(s);
    return s;
}

GLuint BuildProgram(const char* vs, const char* fs)
{
    GLuint v = CompileShader(GL_VERTEX_SHADER,   vs);
    GLuint f = CompileShader(GL_FRAGMENT_SHADER, fs);
    GLuint p = glCreateProgram();
    glAttachShader(p, v);
    glAttachShader(p, f);
    glLinkProgram(p);
    glDeleteShader(v);
    glDeleteShader(f);
    return p;
}

struct WireMesh { vector<glm::vec3> verts; vector<int> lineStarts; vector<int> lineCounts; };

WireMesh GenerateKerrSurface(float r, float spinA, int latLines, int lonLines, int seg)
{
    WireMesh m;
    float rho = sqrtf(r*r + spinA*spinA);
    for (int i = 1; i < latLines; i++) {
        float theta = (float)PI * i / latLines;
        float st = sinf(theta), ct = cosf(theta);
        m.lineStarts.push_back((int)m.verts.size());
        for (int j = 0; j <= seg; j++) {
            float phi = 2.f*(float)PI*j/seg;
            m.verts.push_back({ rho*st*cosf(phi), rho*st*sinf(phi), r*ct });
        }
        m.lineCounts.push_back(seg+1);
    }
    for (int j = 0; j < lonLines; j++) {
        float phi = 2.f*(float)PI*j/lonLines;
        m.lineStarts.push_back((int)m.verts.size());
        for (int i = 0; i <= seg; i++) {
            float theta = (float)PI*i/seg;
            float st = sinf(theta), ct = cosf(theta);
            m.verts.push_back({ rho*st*cosf(phi), rho*st*sinf(phi), r*ct });
        }
        m.lineCounts.push_back(seg+1);
    }
    return m;
}

WireMesh GenerateErgosphere(float M, float spinA, int latLines, int lonLines, int seg)
{
    WireMesh m;
    auto ergR = [&](float theta) {
        float c = cosf(theta);
        float disc = M*M - spinA*spinA*c*c;
        return (disc < 0.f) ? M : M + sqrtf(disc);
    };
    for (int i = 1; i < latLines; i++) {
        float theta = (float)PI * i / latLines;
        float r   = ergR(theta);
        float rho = sqrtf(r*r + spinA*spinA);
        float st  = sinf(theta), ct = cosf(theta);
        m.lineStarts.push_back((int)m.verts.size());
        for (int j = 0; j <= seg; j++) {
            float phi = 2.f*(float)PI*j/seg;
            m.verts.push_back({ rho*st*cosf(phi), rho*st*sinf(phi), r*ct });
        }
        m.lineCounts.push_back(seg+1);
    }
    for (int j = 0; j < lonLines; j++) {
        float phi = 2.f*(float)PI*j/lonLines;
        m.lineStarts.push_back((int)m.verts.size());
        for (int i = 0; i <= seg; i++) {
            float theta = (float)PI*i/seg;
            float r   = ergR(theta);
            float rho = sqrtf(r*r + spinA*spinA);
            float st  = sinf(theta), ct = cosf(theta);
            m.verts.push_back({ rho*st*cosf(phi), rho*st*sinf(phi), r*ct });
        }
        m.lineCounts.push_back(seg+1);
    }
    return m;
}

vector<glm::vec3> GenerateRingSingularity(float spinA, int seg)
{
    vector<glm::vec3> v;
    for (int i = 0; i <= seg; i++) {
        float phi = 2.f*(float)PI*i/seg;
        v.push_back({ spinA*cosf(phi), spinA*sinf(phi), 0.f });
    }
    return v;
}

void DrawWireMesh(const WireMesh& m)
{
    for (int k = 0; k < (int)m.lineStarts.size(); k++)
        glDrawArrays(GL_LINE_STRIP, m.lineStarts[k], m.lineCounts[k]);
}

void MakeVAOVBO(GLuint& vao, GLuint& vbo,
                const void* data = nullptr, size_t bytes = 0,
                GLenum usage = GL_STATIC_DRAW)
{
    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    if (data && bytes)
        glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)bytes, data, usage);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), nullptr);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);
}

static void DataRange(const vector<float>& v, float& lo, float& hi)
{
    lo =  numeric_limits<float>::max();
    hi = -numeric_limits<float>::max();
    for (float x : v) { lo = min(lo, x); hi = max(hi, x); }
    float pad = (hi - lo) == 0.f ? 0.1f : (hi - lo) * 0.10f;
    lo -= pad;  hi += pad;
}

static bool WorldToScreen(const glm::vec3& worldPos,
                           const glm::mat4& MVP,
                           int fbW, int fbH,
                           ImVec2& out)
{
    glm::vec4 clip = MVP * glm::vec4(worldPos, 1.f);
    if (clip.w <= 0.f) return false;
    glm::vec3 ndc = glm::vec3(clip) / clip.w;
    out.x = (ndc.x * 0.5f + 0.5f) * (float)fbW;
    out.y = (1.f - (ndc.y * 0.5f + 0.5f)) * (float)fbH;
    return (ndc.x >= -1.f && ndc.x <= 1.f &&
            ndc.y >= -1.f && ndc.y <= 1.f);
}

int main()
{
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_SAMPLES, 4);

    GLFWwindow* window = glfwCreateWindow(1900, 1000, "GEOstar - RK2", nullptr, nullptr);
    glfwMakeContextCurrent(window);
    gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window,   scroll_callback);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_MULTISAMPLE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");

    ImGuiStyle& style = ImGui::GetStyle();
    style.WindowRounding   = 6.f;
    style.FrameRounding    = 4.f;
    style.GrabRounding     = 4.f;
    style.WindowBorderSize = 0.f;
    style.FramePadding     = {8.f, 5.f};
    ImVec4* c = style.Colors;
    c[ImGuiCol_WindowBg]        = {0.08f, 0.08f, 0.10f, 0.92f};
    c[ImGuiCol_FrameBg]         = {0.14f, 0.14f, 0.18f, 1.f};
    c[ImGuiCol_SliderGrab]      = {0.36f, 0.60f, 1.00f, 1.f};
    c[ImGuiCol_SliderGrabActive]= {0.55f, 0.75f, 1.00f, 1.f};
    c[ImGuiCol_Button]          = {0.20f, 0.35f, 0.65f, 1.f};
    c[ImGuiCol_ButtonHovered]   = {0.30f, 0.50f, 0.85f, 1.f};
    c[ImGuiCol_ButtonActive]    = {0.15f, 0.25f, 0.55f, 1.f};
    c[ImGuiCol_TitleBgActive]   = {0.12f, 0.20f, 0.40f, 1.f};
    c[ImGuiCol_CheckMark]       = {0.36f, 0.60f, 1.00f, 1.f};


    double M = 1.0;
    double a = 0.5;
    KerrMetric Kerr(a, M);
    double Rs = 2 * M;

    vector<StateVector> Bodies(5);

    Bodies[0].X[0]=0;
    Bodies[0].X[1]=6.2;
    Bodies[0].X[2]=PI/2;
    Bodies[0].X[3]=0;
    Bodies[0].X[4]=1.08;
    Bodies[0].X[5]=0.005;
    Bodies[0].X[6]=0.020;
    Bodies[0].X[7]=0.05;

    Bodies[1].X[0]=0;
    Bodies[1].X[1]=6.22;
    Bodies[1].X[2]=PI/3;
    Bodies[1].X[3]=0.02;
    Bodies[1].X[4]=1.08;
    Bodies[1].X[5]=0.006;
    Bodies[1].X[6]=0.021;
    Bodies[1].X[7]=0.05;

    Bodies[2].X[0]=0;
    Bodies[2].X[1]=6.24;
    Bodies[2].X[2]=PI/4;
    Bodies[2].X[3]=0.04;
    Bodies[2].X[4]=1.08;
    Bodies[2].X[5]=0.004;
    Bodies[2].X[6]=0.019;
    Bodies[2].X[7]=0.05;

    Bodies[3].X[0]=0;
    Bodies[3].X[1]=6.26;
    Bodies[3].X[2]=PI/5;
    Bodies[3].X[3]=0.06;
    Bodies[3].X[4]=1.08;
    Bodies[3].X[5]=-0.005;
    Bodies[3].X[6]=0.022;
    Bodies[3].X[7]=0.05;

    Bodies[4].X[0]=0;
    Bodies[4].X[1]=6.28;
    Bodies[4].X[2]=PI/6;
    Bodies[4].X[3]=0.08;
    Bodies[4].X[4]=1.08;
    Bodies[4].X[5]=0.003;
    Bodies[4].X[6]=0.018;
    Bodies[4].X[7]=0.089;

    vector<SolverRES> Systems;
    for (auto& b : Bodies)
        Systems.push_back(RK2Solver(Kerr, b, 1e-3, 200000));

    vector<vector<glm::vec3>> Orbits(5);
    for (int i = 0; i < 5; i++)
        for (auto& X : Systems[i].x)
            Orbits[i].push_back(ToCartesian(X, a));

    int N0 = Systems[0].n;

    vector<float> tau(N0), dtau(N0), vr(N0), vphi(N0);
    for (int i = 0; i < N0; i++) {
        tau[i]  = (float)Systems[0].tau[i];
        dtau[i] = (float)Systems[0].x[i].X[4];
        vr[i]   = (float)Systems[0].x[i].X[5];
        vphi[i] = (float)Systems[0].x[i].X[7];
    }

    float dtauLo, dtauHi, vrLo, vrHi, vphiLo, vphiHi;
    DataRange(dtau, dtauLo, dtauHi);
    DataRange(vr,   vrLo,   vrHi);
    DataRange(vphi, vphiLo, vphiHi);


    const char* VS = R"(
#version 330 core
layout(location=0) in vec3 aPos;
uniform mat4 MVP;
void main() { gl_Position = MVP * vec4(aPos, 1.0); }
)";
    const char* FS = R"(
#version 330 core
out vec4 FragColor;
uniform vec4 uColor;
void main() { FragColor = uColor; }
)";

    GLuint prog     = BuildProgram(VS, FS);
    GLuint locMVP   = glGetUniformLocation(prog, "MVP");
    GLuint locColor = glGetUniformLocation(prog, "uColor");


    GLuint OrbitVAO[5], OrbitVBO[5];
    for (int i = 0; i < 5; i++)
        MakeVAOVBO(OrbitVAO[i], OrbitVBO[i],
                   Orbits[i].data(),
                   Orbits[i].size() * sizeof(glm::vec3),
                   GL_STATIC_DRAW);

    GLuint DotVAO, DotVBO;
    MakeVAOVBO(DotVAO, DotVBO, nullptr, sizeof(glm::vec3), GL_DYNAMIC_DRAW);


    float disc   = (float)sqrt(M*M - a*a);
    float rPlus  = (float)(M + disc);
    float rMinus = (float)(M - disc);

    WireMesh OuterHorizon = GenerateKerrSurface(rPlus,  (float)a, 36, 48, 128);
    WireMesh InnerHorizon = GenerateKerrSurface(rMinus, (float)a, 24, 32, 128);
    WireMesh Ergo         = GenerateErgosphere((float)M, (float)a, 36, 48, 128);
    vector<glm::vec3> RingSing = GenerateRingSingularity((float)a, 256);

    GLuint OuterVAO, OuterVBO;
    MakeVAOVBO(OuterVAO, OuterVBO,
               OuterHorizon.verts.data(),
               OuterHorizon.verts.size() * sizeof(glm::vec3));

    GLuint InnerVAO, InnerVBO;
    MakeVAOVBO(InnerVAO, InnerVBO,
               InnerHorizon.verts.data(),
               InnerHorizon.verts.size() * sizeof(glm::vec3));

    GLuint ErgoVAO, ErgoVBO;
    MakeVAOVBO(ErgoVAO, ErgoVBO,
               Ergo.verts.data(),
               Ergo.verts.size() * sizeof(glm::vec3));

    GLuint RingVAO, RingVBO;
    MakeVAOVBO(RingVAO, RingVBO,
               RingSing.data(),
               RingSing.size() * sizeof(glm::vec3));




    auto ergPoleR = [&]() -> float {

        return (float)(M + sqrt(M*M - a*a));
    };


    glm::vec3 ergoLabelAnchor  = { 0.f, 0.f, ergPoleR() };
    glm::vec3 outerLabelAnchor = { 0.f, 0.f, rPlus       };
    glm::vec3 innerLabelAnchor = { 0.f, 0.f, rMinus       };


    int   currentStep = 0;
    float speed       = 300.f;
    bool  paused      = false;
    bool  bodyVisible[5] = {true,true,true,true,true};

    float prev = (float)glfwGetTime();


    while (!glfwWindowShouldClose(window))
    {
        glfwPollEvents();

        float now = (float)glfwGetTime();
        float dt  = now - prev;
        prev = now;


        if (!paused) {
            currentStep += (int)(speed * dt);
            if (currentStep >= N0) currentStep = N0 - 1;
        }

        int plotCount = max(2, currentStep + 1);

        int fbW, fbH;
        glfwGetFramebufferSize(window, &fbW, &fbH);
        glViewport(0, 0, fbW, fbH);

        glClearColor(0.04f, 0.04f, 0.06f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        mouseActive = !ImGui::GetIO().WantCaptureMouse;


        glm::vec3 cam = {
            radius * cosf(glm::radians(yaw))   * cosf(glm::radians(pitch)),
            radius * sinf(glm::radians(pitch)),
            radius * sinf(glm::radians(yaw))   * cosf(glm::radians(pitch))
        };
        glm::mat4 view = glm::lookAt(cam, {0,0,0}, {0,1,0});
        glm::mat4 proj = glm::perspective(glm::radians(45.f),
                                          (float)fbW / (float)fbH,
                                          0.1f, 1000.f);
        glm::mat4 MVP = proj * view;


        ImGui::SetNextWindowPos({10.f, 10.f}, ImGuiCond_Always);
        ImGui::SetNextWindowSize({310.f, 250.f}, ImGuiCond_Always);
        ImGui::Begin("Controls", nullptr,
            ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
            ImGuiWindowFlags_NoCollapse);
        ImGui::SliderFloat("Speed", &speed, 10.f, 2000.f, "%.0f steps/s");
        ImGui::Spacing();
        if (paused) { if (ImGui::Button("  Resume  ")) paused = false; }
        else        { if (ImGui::Button("  Pause   ")) paused = true;  }
        ImGui::SameLine();
        if (ImGui::Button("  Reset   ")) { currentStep = 0; paused = false; }
        ImGui::Spacing();
        ImGui::Separator();
        ImGui::TextDisabled("Body visibility");
        for (int i = 0; i < 5; i++) {
            ImGui::PushStyleColor(ImGuiCol_CheckMark,
                ImVec4(OrbitColors[i].x, OrbitColors[i].y, OrbitColors[i].z, 1.f));
            string lbl = "Body " + to_string(i + 1);
            ImGui::Checkbox(lbl.c_str(), &bodyVisible[i]);
            ImGui::PopStyleColor();
            if (i < 4) ImGui::SameLine();
        }
        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();
        ImGui::TextDisabled("Kerr  a=%.2f  M=%.2f  Radius of the Ergoshere=%.2f", a, M, Rs);
        ImGui::TextDisabled("r+ =%.3f  r- =%.3f", rPlus, rMinus);
        ImGui::End();


        ImGui::SetNextWindowPos({(float)fbW - 380.f, 10.f}, ImGuiCond_Always);
        ImGui::SetNextWindowSize({370.f, 200.f}, ImGuiCond_Always);
        ImGui::Begin("Coordinates", nullptr,
            ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
            ImGuiWindowFlags_NoCollapse);
        for (int i = 0; i < 5; i++) {
            int idx = min(currentStep, Systems[i].n - 1);
            ImGui::PushStyleColor(ImGuiCol_Text,
                ImVec4(OrbitColors[i].x, OrbitColors[i].y, OrbitColors[i].z, 1.f));
            ImGui::Text("Body %d", i + 1);
            ImGui::PopStyleColor();
            ImGui::SameLine(75.f);
            ImGui::Text("r = %7.3f", Systems[i].x[idx].X[1]);
            ImGui::SameLine(195.f);
            ImGui::Text("phi = %6.3f", Systems[i].x[idx].X[3]);
            if (i < 4) ImGui::Separator();
        }
        ImGui::End();


        ImGui::SetNextWindowPos({10.f, 270.f}, ImGuiCond_Always);
        ImGui::SetNextWindowSize({210.f, 160.f}, ImGuiCond_Always);
        ImGui::Begin("Legend", nullptr,
            ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
            ImGuiWindowFlags_NoCollapse);


        auto LegendRow = [&](ImVec4 col, const char* label) {
            ImDrawList* dl = ImGui::GetWindowDrawList();
            ImVec2 p = ImGui::GetCursorScreenPos();
            dl->AddRectFilled({p.x, p.y + 2.f}, {p.x + 14.f, p.y + 14.f}, ImGui::ColorConvertFloat4ToU32(col), 2.f);
            ImGui::Dummy({16.f, 16.f});
            ImGui::SameLine();
            ImGui::TextUnformatted(label);
        };

        LegendRow({0.95f, 0.75f, 0.10f, 0.80f}, "Ergosphere");
        LegendRow({0.90f, 0.25f, 0.05f, 0.90f}, "Outer horizon (r+)");
        LegendRow({0.35f, 0.55f, 1.00f, 0.90f}, "Inner horizon (r-)");
        LegendRow({1.00f, 1.00f, 1.00f, 1.00f}, "Ring singularity");
        ImGui::End();


        float diagH = 220.f;
        float diagW = (float)fbW - 20.f;
        ImGui::SetNextWindowPos({10.f, (float)fbH - diagH - 10.f}, ImGuiCond_Always);
        ImGui::SetNextWindowSize({diagW, diagH}, ImGuiCond_Always);
        ImGui::Begin("Diagnostics (Blue Body)", nullptr,
            ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
            ImGuiWindowFlags_NoCollapse);

        float xMin = tau[0];
        float xMax = tau[N0 - 1];

        ImVec2 plotSz = { (diagW - 40.f) / 3.f, 155.f };

        if (ImPlot::BeginPlot("Local Time Dilation -dt/dtau- ", plotSz)) {
            ImPlot::SetupAxes("tau", "u_t");
            ImPlot::SetupAxisLimits(ImAxis_X1, xMin, xMax, ImGuiCond_Always);
            ImPlot::SetupAxisLimits(ImAxis_Y1, dtauLo, dtauHi, ImGuiCond_Always);
            ImPlot::PlotLine("u_t", tau.data(), dtau.data(), plotCount);
            ImPlot::EndPlot();
        }
        ImGui::SameLine();

        if (ImPlot::BeginPlot("Radial speed ", plotSz)) {
            ImPlot::SetupAxes("tau", "u_r");
            ImPlot::SetupAxisLimits(ImAxis_X1, xMin, xMax, ImGuiCond_Always);
            ImPlot::SetupAxisLimits(ImAxis_Y1, vrLo, vrHi, ImGuiCond_Always);
            ImPlot::PlotLine("u_r", tau.data(), vr.data(), plotCount);
            ImPlot::EndPlot();
        }
        ImGui::SameLine();

        if (ImPlot::BeginPlot("Azimuthal speed ", plotSz)) {
            ImPlot::SetupAxes("tau", "u_phi");
            ImPlot::SetupAxisLimits(ImAxis_X1, xMin, xMax, ImGuiCond_Always);
            ImPlot::SetupAxisLimits(ImAxis_Y1, vphiLo, vphiHi, ImGuiCond_Always);
            ImPlot::PlotLine("u_phi", tau.data(), vphi.data(), plotCount);
            ImPlot::EndPlot();
        }

        ImGui::End();


        glUseProgram(prog);
        glUniformMatrix4fv(locMVP, 1, GL_FALSE, glm::value_ptr(MVP));


        glBindVertexArray(ErgoVAO);
        glUniform4f(locColor, 0.95f, 0.75f, 0.10f, 0.30f);
        DrawWireMesh(Ergo);


        glBindVertexArray(OuterVAO);
        glUniform4f(locColor, 0.90f, 0.25f, 0.05f, 0.70f);
        DrawWireMesh(OuterHorizon);


        glBindVertexArray(InnerVAO);
        glUniform4f(locColor, 0.35f, 0.55f, 1.00f, 0.85f);
        DrawWireMesh(InnerHorizon);


        glBindVertexArray(RingVAO);
        glUniform4f(locColor, 1.00f, 1.00f, 1.00f, 1.00f);
        glDrawArrays(GL_LINE_STRIP, 0, (GLsizei)RingSing.size());


        float bodyPulse = 0.5f + 0.5f * sinf(now * 3.f);
        for (int i = 0; i < 5; i++) {
            if (!bodyVisible[i]) continue;
            glm::vec3 col = OrbitColors[i];
            int count = min(currentStep + 1, (int)Orbits[i].size());

            glBindVertexArray(OrbitVAO[i]);
            glUniform4f(locColor, col.x, col.y, col.z, 0.80f);
            glDrawArrays(GL_LINE_STRIP, 0, count);

            if (currentStep < (int)Orbits[i].size()) {
                glBindVertexArray(DotVAO);
                glBindBuffer(GL_ARRAY_BUFFER, DotVBO);
                glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3),
                             &Orbits[i][currentStep], GL_DYNAMIC_DRAW);
                glUniform4f(locColor, col.x, col.y, col.z, 1.f);
                glPointSize(8.f + 6.f * bodyPulse);
                glDrawArrays(GL_POINTS, 0, 1);
            }
        }


        {
            ImDrawList* bg = ImGui::GetBackgroundDrawList();


            auto DrawLabel3D = [&](const glm::vec3& anchor,
                                   const char* text,
                                   ImVec4 ,
                                   ImVec2 offset = {0.f, -18.f})
            {
                ImVec2 sp;
                if (!WorldToScreen(anchor, MVP, fbW, fbH, sp)) return;
                sp.x += offset.x;
                sp.y += offset.y;

                ImVec2 textSz = ImGui::CalcTextSize(text);
                float pad = 4.f;
                ImU32 bgCol   = IM_COL32(10, 10, 14, 180);
                ImU32 textCol = IM_COL32(210, 210, 210, 230);

                ImVec2 pMin = { sp.x - pad, sp.y - pad };
                ImVec2 pMax = { sp.x + textSz.x + pad, sp.y + textSz.y + pad };
                bg->AddRectFilled(pMin, pMax, bgCol, 3.f);
                bg->AddText({ sp.x, sp.y }, textCol, text);
            };


            DrawLabel3D(ergoLabelAnchor,
                        "Ergosphere",
                        {0.95f, 0.75f, 0.10f, 1.f},
                        {0.f, -22.f});


            DrawLabel3D(outerLabelAnchor,
                        "Outer horizon (r+)",
                        {0.90f, 0.25f, 0.05f, 1.f},
                        {0.f, -22.f});


            DrawLabel3D(innerLabelAnchor,
                        "Inner horizon (r-)",
                        {0.35f, 0.55f, 1.00f, 1.f},
                        {0.f, -18.f});

        }

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }

    ImPlot::DestroyContext();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}