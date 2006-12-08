/*
 * Copyright (C) 2004-2007 Andrew Mihal
 *
 * This file is part of Enblend.
 *
 * Enblend is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * Enblend is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Enblend; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <iostream>

#include "gpu.h"

using std::cerr;
using std::cout;
using std::endl;

static GLuint GlutWindowHandle;
static int MaxTextureSize;

void checkGLErrors(int line, char *file) {
    GLenum errCode;
    if ((errCode = glGetError()) != GL_NO_ERROR) {
        cerr << "enblend: GL error in " << file << ":" << line << ": " << gluErrorString(errCode) << endl;
        exit(1);
    }
}

void printInfoLog(GLhandleARB obj) {
    int infologLength = 0;
    int charsWritten = 0;
    char *infoLog;
    glGetObjectParameterivARB(obj, GL_OBJECT_INFO_LOG_LENGTH_ARB, &infologLength);
    if (infologLength > 1) {
        infoLog = new char[infologLength];
        glGetInfoLogARB(obj, infologLength, &charsWritten, infoLog);
        cout << "enblend: GL info log:" << endl << infoLog << endl;
        delete[] infoLog;
    }
}

bool checkFramebufferStatus() {
    GLenum status;
    status = (GLenum) glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
    switch(status) {
        case GL_FRAMEBUFFER_COMPLETE_EXT:
            return true;
        case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT_EXT:
            cerr << "enblend: GL error: Framebuffer incomplete, incomplete attachment" << endl;
            return false;
        case GL_FRAMEBUFFER_UNSUPPORTED_EXT:
            cerr << "enblend: Unsupported framebuffer format" << endl;
            return false;
        case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT_EXT:
            cerr << "enblend: Framebuffer incomplete, missing attachment" << endl;
            return false;
        case GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS_EXT:
            cerr << "enblend: Framebuffer incomplete, attached images must have same dimensions" << endl;
            return false;
        case GL_FRAMEBUFFER_INCOMPLETE_FORMATS_EXT:
            cerr << "enblend: Framebuffer incomplete, attached images must have same format" << endl;
            return false;
        case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER_EXT:
            cerr << "enblend: Framebuffer incomplete, missing draw buffer" << endl;
            return false;
        case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER_EXT:
            cerr << "enblend: Framebuffer incomplete, missing read buffer" << endl;
            return false;
    }

    return false;
}

bool initGPU() {
    glutInitDisplayMode(GLUT_RGBA | GLUT_ALPHA);
    GlutWindowHandle = glutCreateWindow("Enblend");

    int err = glewInit();
    if (err != GLEW_OK) {
        cerr << "enblend: an error occured while setting up the GPU:" << endl;
        cerr << glewGetErrorString(err) << endl;
        cerr << "enblend: sorry, the --gpu flag is not going to work on this machine." << endl;
        glutDestroyWindow(GlutWindowHandle);
        exit(1);
    }

    cout << "enblend: using graphics card: " << glGetString(GL_VENDOR) << " " << glGetString(GL_RENDERER) << endl;

    bool has_arb_fragment_shader = glewGetExtension("GL_ARB_fragment_shader");
    bool has_arb_vertex_shader = glewGetExtension("GL_ARB_vertex_shader");
    bool has_arb_shader_objects = glewGetExtension("GL_ARB_shader_objects");
    bool has_arb_shading_language = glewGetExtension("GL_ARB_shading_language_100");

    cout << "enblend: checking extensions: GL_ARB_fragment_shader = " << has_arb_fragment_shader << endl;
    cout << "enblend: checking extensions: GL_ARB_vertex_shader = " << has_arb_vertex_shader << endl;
    cout << "enblend: checking extensions: GL_ARB_shader_objects = " << has_arb_shader_objects << endl;
    cout << "enblend: checking extensions: GL_ARB_shading_language_100 = " << has_arb_shading_language << endl;

    if (!(has_arb_fragment_shader && has_arb_vertex_shader && has_arb_shader_objects && has_arb_shading_language)) {
        cerr << "enblend: this graphics card lacks the necessary extensions for --gpu." << endl;
        cerr << "enblend: sorry, the --gpu flag is not going to work on this machine." << endl;
        glutDestroyWindow(GlutWindowHandle);
        exit(1);
    }

    glGetIntegerv(GL_MAX_TEXTURE_SIZE, &MaxTextureSize);

    return true;
}

bool wrapupGPU() {
    glutDestroyWindow(GlutWindowHandle);
    return true;
}

