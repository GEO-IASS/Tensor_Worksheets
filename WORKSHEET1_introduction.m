%% Introduction to the worksheets
% _Written by Sebastian Kraemer, IGPM at RWTH Aachen University_
%
% You can *view* this document properly using the Publish tool, which will open
% a .html file (RUN_ME_FIRST.m may just have done this for you). You can find 
% that option under PUBLISH (right next to the EDITOR tab above). Doing so will also
% display a table of contents. Best try it out now.
%
% The worksheets will be Matlab files with inline code. We will not
% test your knowledge about the *exercises* (we surely can if you
% wish), but nonetheless it is a good idea to treat them as if we would.
% There will be small sections in which you can enter your *solutions*. You
% can of course modify anything you like, but in that case make sure that
% the exercises remain intact.
%
% This worksheet contains some *help* for Matlab, but if you have trouble or
% just never worked with Matlab, just ask for a short introduction. At any
% time, it is always a good idea to look into the very good documentation
% of Matlab. If this does not work out, you may just write Pseudo Code and ask for a
% translation into Matlab.
%
% To *evaluate* a single section (in between horizontal lines) click any line
% within the section and press alt+enter. You may also do this via a right
% click. The result will appear in the command line window.
%
% You can use the command line window also to *play around*, for example to
% understand a new function. For longer tests, or tests you like to keep for a while,
% use SANDBOX.m. 
%
% It may further be helpful to change the default *shortcuts* from emacs to
% windows (select the HOME tab above and navigate to
% Preferences/Keyboard/Shortcuts).
%%
clear all % this clears all variables of their values

%% Constructing a tensor
% Luckily, Matlab already has tensors at its disposel, which work very
% analogous to matrices. Not suprisingly, it is slightly more complicated to
% display these.
A1 = zeros(2,3); % a semicolon suppresses the output
A1 % only now is A1 displayed
T1 = zeros(2,3,4)

%% Implemention of tensors
% Matlab (and any other code), at its most basic level, saves tensors as array (or vectors). 
% So, every time the operator (i_1,...,i_d) is called, the computer needs to find the 
% respective entry in this array. Matlab can turn any tensor to this array
% using the (:) command. This is very helpful to fill a tensor. As you
% will see, matrices are saved column-wise, while for tensor mode 1
% comes first, too, and is then followed by mode 2, 3, ...:
A2 = [1,3,5; 2,4,6]
A2(:)
%%
T1(:) = 1:(2*3*4)
T1(1,2,3)
T1(15)

%% Using a for loop
% One may also use for loops to fill in the entries of a tensor:
B1 = zeros(3,4,2);
for i = 1:3*4*2
    B1(i) = sin(pi/2*i);
end
B1(:,:,1)

% or simply use B1(:) = sin(1:4*3*2)

%% EXERCISE 1: operator()
% Construct a 3 by 2 by 4 tensor (default value 0) and set entry number 12 as 1.
% To which index (i_1,...,i_d) of the tensor does this entry belong?
% Call the corresponding Matlab operator to verify your answer.

%% SOLUTION 1
T2 = 'ANSWER 2 MISSING';

%% The reshape function
% Tensors can be reshaped without any relevant cost, since Matlab does
% not need to change the underlying array of a tensor. Surely, this
% operation must not change the total number of entries.
size(T1)
T1(15)
T1(2,2,2)
%%
T3 = reshape(T1,[3,4,2]);
size(T3)
T3(15)
T3(2,2,2)

% reshape(T1,[1,2,3]) produces an error

%% EXERCISE 2: reshape
% Explain what happens in the following line:
A2 = reshape(T1,[6,4])

%% EXERCISE 3: operator_bracket
% Complete the operator_bracket function which should do exactly the same as the operator
% (i_1,...,i_d), that is, return the entry T(i_1,...,i_d), where I = [i_1,...,i_d]
% (but of course without using that operator).
% You may select the function below with the cursor and press control+d to 
% open the separate file operator_bracket directly. 
T4 = rand(2,3,4,5);
T4(1,2,3,4)
operator_bracket(T4,[1,2,3,4])
T4(2,1,2,1)
operator_bracket(T4,[2,1,2,1])

%% Next: WORKSHEET2_CP_format.m
% open WORKSHEET2_CP_format.m and continue there













