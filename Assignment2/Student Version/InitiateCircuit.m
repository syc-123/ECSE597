% Author: Karan 
% Date: 13/09/2022

% This script initiates the structures and map containers for all the
% element types.

clear all
global elementList nodeMap

% elements Supported


nodeMap = containers.Map('KeyType','char','ValueType','int32');

elementList.Resistors = {};
elementList.Resistors.containerMap = containers.Map('KeyType','char','ValueType','int32');
elementList.Resistors.numElements =0 ;

elementList.Capacitors = {};
elementList.Capacitors.containerMap = containers.Map('KeyType','char','ValueType','int32');
elementList.Capacitors.numElements =0 ;

elementList.Inductors = {};
elementList.Inductors.containerMap = containers.Map('KeyType','char','ValueType','int32');
elementList.Inductors.numElements =0 ;


% voltage controlled current sources
elementList.VCCS ={};
elementList.VCCS.containerMap = containers.Map('KeyType','char','ValueType','int32');
elementList.VCCS.numElements =0 ;

% voltage controlled voltage sources
elementList.VCVS ={};
elementList.VCVS.containerMap = containers.Map('KeyType','char','ValueType','int32');
elementList.VCVS.numElements =0 ;


% current controlled current sources
elementList.CCCS ={};
elementList.CCCS.containerMap = containers.Map('KeyType','char','ValueType','int32');
elementList.CCCS.numElements =0 ;

% current controlled voltsge sources
elementList.CCVS ={};
elementList.CCVS.containerMap = containers.Map('KeyType','char','ValueType','int32');
elementList.CCVS.numElements =0 ;

% add voltage sources
elementList.DC_VolSources = {};
elementList.DC_VolSources.containerMap = containers.Map('KeyType','char','ValueType','int32');
elementList.DC_VolSources.numElements =0 ;

elementList.AC_VolSources = {};
elementList.AC_VolSources.containerMap = containers.Map('KeyType','char','ValueType','int32');
elementList.AC_VolSources.numElements =0 ;

elementList.TR_VolSources = {};
elementList.TR_VolSources.containerMap = containers.Map('KeyType','char','ValueType','int32');
elementList.TR_VolSources.numElements =0 ;

elementList.HB_VolSources = {};
elementList.HB_VolSources.containerMap = containers.Map('KeyType','char','ValueType','int32');
elementList.HB_VolSources.numElements =0 ;
 
% add current sources
elementList.DC_CurSources = {};
elementList.DC_CurSources.containerMap = containers.Map('KeyType','char','ValueType','int32');
elementList.DC_CurSources.numElements =0 ;

elementList.AC_CurSources = {};
elementList.AC_CurSources.containerMap = containers.Map('KeyType','char','ValueType','int32');
elementList.AC_CurSources.numElements =0 ;


elementList.TR_CurSources = {};
elementList.TR_CurSources.containerMap = containers.Map('KeyType','char','ValueType','int32');
elementList.TR_CurSources.numElements =0 ;


elementList.HB_CurSources = {};
elementList.HB_CurSources.containerMap = containers.Map('KeyType','char','ValueType','int32');
elementList.HB_CurSources.numElements =0 ;

% Opamp
elementList.OpAmps = {};
elementList.OpAmps.containerMap = containers.Map('KeyType','char','ValueType','int32');
elementList.OpAmps.numElements =0 ;

% Mutual Inductors 
elementList.Mutual_Inductors = {};
elementList.Mutual_Inductors.containerMap = containers.Map('KeyType','char','ValueType','int32');
elementList.Mutual_Inductors.numElements =0 ;


% Diodes
elementList.Diodes = {};
elementList.Diodes.containerMap = containers.Map('KeyType','char','ValueType','int32');
elementList.Diodes.numElements =0 ;

% BJTs
elementList.BJTs = {};
elementList.BJTs.containerMap = containers.Map('KeyType','char','ValueType','int32');
elementList.BJTs.numElements =0 ;

